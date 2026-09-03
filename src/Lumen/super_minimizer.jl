# ============================================================================ #
#   SUPER LUMEN: seal-then-combine, single-pass rule extraction                #
#                                                                              #
#
#   PATCH (questa revisione): STRATEGIA "SEAL-THEN-COMBINE" AL POSTO DEL      #
#   FOLD INCREMENTALE.                                                        #
#   ------------------------------------------------------------------
#   La revisione precedente faceva, per ogni classe, un fold INCREMENTALE:
#
#       φ  + raw_batch₁ → φ₁
#       φ₁ + raw_batch₂ → φ₂
#       φ₂ + raw_batch₃ → φ₃
#       ...
#
#   cioè ogni volta che il buffer raggiungeva M righe, Espresso veniva
#   invocato su (termini-già-minimizzati-finora ∪ nuovo raw), e il
#   risultato RIMPIAZZAVA la formula corrente. La formula "corrente" era
#   quindi sempre una sola, che cresceva/si comprimeva ad ogni fold.
#
#   Questa revisione cambia la strategia in "seal-then-combine":
#
#       raw_batch₁ → φ₁   (sigillato, MESSO DA PARTE, mai più toccato)
#       raw_batch₂ → φ₂   (sigillato, MESSO DA PARTE, mai più toccato)
#       raw_batch₃ → φ₃   (sigillato, MESSO DA PARTE, mai più toccato)
#       ...
#       (fine pass) → combinazione: Φ = φ₁ ∪ φ₂ ∪ φ₃ ∪ ...
#       (opzionale) → una SOLA ri-minimizzazione finale di Φ
#
#   Cioè: ogni batch di raw cube viene minimizzato ISOLATAMENTE (mai contro
#   una φ precedente), il risultato di quella minimizzazione è una φᵢ
#   indipendente che viene semplicemente accantonata in una lista
#   (`buf.phi_groups`) e MAI più data in pasto a Espresso insieme al
#   prossimo batch. Solo a pass conclusa (tutte le classi, tutto lo spazio
#   enumerato o tutti i fold intermedi scattati) si fa la "combine": si
#   concatenano tutti i φᵢ di ciascuna classe (`vcat`, che essendo pura
#   unione di ON-set-coprenti non può mai perdere copertura) e, se
#   richiesto (`reminimize_after_combine=true`, default), si fa UNA sola
#   passata finale di minimizzazione sul totale per rimuovere la
#   ridondanza lasciata dal sigillo "a pezzi".
#
#   Perché farlo: nel fold incrementale, ogni volta che Espresso rivede la
#   φ corrente insieme a nuovo raw, può "rimescolare" termini che erano
#   già stabili, cambiando la forma della formula ad ogni fold in modi che
#   dipendono fortemente dall'ordine di enumerazione random. Con
#   seal-then-combine, ogni φᵢ è definitiva nel momento in cui viene
#   creata: non viene mai ri-derivata a partire da una φ precedente, quindi
#   il contributo di ciascun batch alla formula finale è più isolato e
#   riproducibile batch-per-batch (a parità di offset/OFF-set), al prezzo
#   di rimandare tutta la deduplicazione fra batch diversi alla combine
#   finale.
#
#   NOTA IMPORTANTE SUL BUDGET M: con il fold incrementale, il controllo
#   "ok = length(terms) < M dopo il fold" era una garanzia (euristica) che
#   la formula CUMULATIVA restasse sotto M. Con seal-then-combine questo
#   NON è più vero: il controllo `< M` ora si applica SOLO alla singola φᵢ
#   appena sigillata (cioè: "questo batch, minimizzato da solo, produce
#   una φᵢ ragionevolmente piccola"), non alla somma di tutte le φᵢ già
#   accantonate. È una conseguenza intrinseca dello "sigilla ora, combina
#   dopo": se serve un vero e proprio hard-cap sul totale, va applicato
#   DOPO la combine finale (fuori dal loop di enumerazione), non durante.
#
#   Tutto il resto (OFF-set esplicito via `.type fr` per evitare che
#   Espresso tratti come OFF territorio non ancora enumerato o di altre
#   classi, LazyPermutation, derivazione soglie/feature, loop di restart
#   con scoring) resta come nella revisione precedente -- vedi la vecchia
#   nota di patch qui sotto per il dettaglio di quel fix, che è ancora
#   quello che rende sound l'OFF-set passato a ogni sigillo/combine.
#
# ---------------------------------------------------------------------------- #
#   (nota di patch precedente, ancora valida per la logica dell'OFF-set)
#
#   Root cause identificata a suo tempo: passare a Espresso solo l'ON-set
#   (PLA `.type f` implicito) lo induceva a trattare come OFF ogni riga non
#   ancora enumerata per QUELLA classe -- incluse righe della stessa
#   classe non ancora viste in questo restart, e territorio di ALTRE
#   classi. Fix: OFF-set esplicito (`.type fr`, vedi `PLA_patch.jl`),
#   costruito come union dei cubi di tutte le ALTRE classi osservate finora
#   (sound perché il modello è deterministico: una riga -> una sola
#   classe). Tutto quello che non è stato ancora enumerato per NESSUNA
#   classe resta don't-care implicito.
# ============================================================================ #

# ---------------------------------------------------------------------------- #
#                              debug instrumentation                           #
# ---------------------------------------------------------------------------- #
const _SEQ_DEBUG = Ref(false)
LUMEN_SEQ_DEBUG!(flag::Bool) = (_SEQ_DEBUG[]=flag; nothing)

@inline function _dbg(args...)
    if _SEQ_DEBUG[]
        println("[super_lumen] ", args...)
    end
    return nothing
end


# ---------------------------------------------------------------------------- #
#                     LAZY RANDOM PERMUTATION (no materialization)             #
# ---------------------------------------------------------------------------- #
struct LazyPermutation
    n::Int
    half_bits::Int
    mask::UInt64
    rounds::Int
    key::UInt64
end

function LazyPermutation(n::Integer; seed::Integer=rand(UInt64), rounds::Int=4)
    n <= 0 && throw(ArgumentError("n must be positive, got $n"))
    total_bits = max(2, ceil(Int, log2(n)))
    total_bits = iseven(total_bits) ? total_bits : total_bits + 1
    half_bits = total_bits ÷ 2
    mask = (UInt64(1) << half_bits) - UInt64(1)
    return LazyPermutation(Int(n), half_bits, mask, rounds, UInt64(seed))
end

@inline function _feistel_round(l::UInt64, r::UInt64, round::Int, p::LazyPermutation)
    h = (hash((r, round, p.key)) & p.mask)
    newl = r
    newr = (l ⊻ h) & p.mask
    return newl, newr
end

@inline function _feistel_permute(x::UInt64, p::LazyPermutation)
    l = (x >> p.half_bits) & p.mask
    r = x & p.mask
    for round in 1:p.rounds
        l, r = _feistel_round(l, r, round, p)
    end
    return (l << p.half_bits) | r
end

function permute(p::LazyPermutation, i::Integer)
    (0 <= i < p.n) || throw(BoundsError(p, i))
    x = UInt64(i)
    while true
        x = _feistel_permute(x, p)
        x < UInt64(p.n) && return Int(x)
    end
end


# ---------------------------------------------------------------------------- #
#                    O(1) linear-index -> combo tuple decoding                 #
# ---------------------------------------------------------------------------- #
function _strides(lens::Vector{Int})
    n = length(lens)
    strides = ones(Int, n)
    @inbounds for j in 2:n
        strides[j] = strides[j-1] * lens[j-1]
    end
    return strides
end

@inline function _combo_at(
    thrs_with_p::Vector{<:AbstractVector},
    lens::Vector{Int},
    strides::Vector{Int},
    linidx0::Int
)
    return ntuple(length(thrs_with_p)) do j
        idx = (linidx0 ÷ strides[j]) % lens[j] + 1
        @inbounds thrs_with_p[j][idx]
    end
end


# ---------------------------------------------------------------------------- #
#                       shared threshold/family derivation                     #
# ---------------------------------------------------------------------------- #
function _prepare_sequential_context(
    config::LumenConfig{T},
    model::SM.AbstractModel
) where {T<:AbstractFloat}
    depth = config.depth

    atoms = unique!(_normalize_atom.(if depth < 1.0
        mapreduce(vcat, SM.models(model); init=SL.Atom{SD.ScalarCondition}[]) do t
            _take_first_percentage(_extract_atoms_bfs_order(t), depth)
        end
    else
        SL.atoms(SM.alphabet(model, false))
    end))

    let unsupported = unique(op for op in get_operator.(atoms) if op ∉ _supported_operators)
        isempty(unsupported) || throw(ArgumentError(
            "Only '<', '≥', '>', '≤' operators are currently supported. " *
            "Found unsupported operators: $(unsupported)."
        ))
    end

    features = SM.featurename.(unique!(get_feature.(atoms)))
    featurenames = SM.info(model, :featurenames)
    classnames = unique!(SM.info(model, :supporting_labels))

    thresholds = Vector{Vector{T}}(undef, length(featurenames))
    op_families = Vector{Symbol}(undef, length(featurenames))

    @inbounds for i in eachindex(featurenames)
        idx = findfirst(f -> f == featurenames[i], features)
        if isnothing(idx)
            thresholds[i] = T[]
            op_families[i] = :lt
        else
            family = _feature_op_family(atoms, features[idx])
            op_families[i] = family
            thresholds[i] = sort!(
                get_threshold.(_atoms_for_feature(atoms, features[idx]));
                rev=(family === :lt)
            )
        end
    end

    thrs_with_p = _thrs_with_boundary(thresholds, op_families)
    lens = length.(thrs_with_p)
    strides = _strides(lens)
    n_total = prod(lens)

    _dbg("context prepared: n_total=", n_total, " nfeatures=", length(featurenames),
        " nclasses=", length(classnames), " lens=", lens)

    return (;
        thresholds, thrs_with_p, featurenames=Symbol.(featurenames),
        classnames, op_families, lens, strides, n_total
    )
end


# ---------------------------------------------------------------------------- #
#                             per-class raw buffer                             #
# ---------------------------------------------------------------------------- #
"""
    _SeqBuffer

Stato per-classe durante un restart.

- `raw`: cubi grezzi ancora NON sigillati, in attesa che il batch
  raggiunga M righe (o che finisca lo spazio da enumerare).
- `phi_groups`: lista di gruppi di termini GIÀ MINIMIZZATI e SIGILLATI.
  Ogni elemento è il risultato di UNA sola chiamata a Espresso su UN solo
  batch di raw cube, presi ISOLATAMENTE. Una volta che un gruppo entra in
  questa lista non viene MAI più ri-minimizzato insieme a un batch
  successivo -- questo è esattamente ciò che distingue "seal-then-combine"
  dal vecchio fold incrementale (che invece rimpiazzava una singola
  formula corrente ad ogni fold).
- `n_rows_seen`: solo per debug/diagnostica, conteggio totale di righe
  osservate per questa classe in questo restart.
"""
mutable struct _SeqBuffer
    raw::Vector{Vector{SL.Atom}}
    phi_groups::Vector{Vector{Any}}
    n_rows_seen::Int
end

_SeqBuffer() = _SeqBuffer(Vector{Vector{SL.Atom}}(), Vector{Vector{Any}}(), 0)

# Occupazione ai fini del trigger di sigillo: SOLO il raw non ancora
# sigillato conta verso M. I φ già sigillati sono "fuori budget" perché
# non vengono più ri-processati da Espresso finché non si arriva alla
# combine finale (fuori dal loop di enumerazione).
@inline _raw_occupancy(buf::_SeqBuffer) = length(buf.raw)

# Solo per debug: quanti termini totali sono già stati sigillati finora
# per questa classe, sommando tutti i gruppi.
@inline _n_sealed_terms(buf::_SeqBuffer) = sum(length, buf.phi_groups; init=0)


# ---------------------------------------------------------------------------- #
#                     coverage / conversion helper utilities                   #
# ---------------------------------------------------------------------------- #
function _cube_satisfies_term(cube::Vector{SL.Atom}, term_atoms::Vector{SL.Atom})
    for a in term_atoms
        cond_a = SL.value(a)
        ok = any(cube) do b
            cond_b = SL.value(b)
            SD.feature(cond_b) == SD.feature(cond_a) &&
                (cond_b == cond_a || SD.includes(cond_b, cond_a))
        end
        ok || return false
    end
    return true
end

function _verify_fold_coverage(
    combined::Vector{Vector{SL.Atom}},
    minimized_cubes::Vector{Vector{SL.Atom}};
    class_label,
    M::Int
)
    _SEQ_DEBUG[] || return true
    uncovered = 0
    for (i, row) in enumerate(combined)
        covered = any(t -> _cube_satisfies_term(row, t), minimized_cubes)
        if !covered
            uncovered += 1
            uncovered <= 3 && _dbg("  COVERAGE LOSS class=", class_label,
                " combined[", i, "]=", row,
                " NOT covered by any minimized term")
        end
    end
    if uncovered > 0
        _dbg("*** COVERAGE CHECK FAILED *** class=", class_label,
            " uncovered=", uncovered, "/", length(combined), " M=", M)
        return false
    else
        _dbg("coverage check OK: class=", class_label,
            " all ", length(combined), " rows covered")
        return true
    end
end

_lit_to_atom(child::SL.Atom) = child
function _lit_to_atom(child::SL.Literal)
    if SL.ispos(child)
        return SL.atom(child)
    end
    cond = SL.value(SL.atom(child))
    if !SD.hasdual(cond)
        throw(ArgumentError(
            "_term_to_cube: negative literal on condition $(cond) has no " *
            "dual condition to rewrite it to -- cannot round-trip this " *
            "term back into unsigned-atom cube form."
        ))
    end
    return SL.Atom(SD.dual(cond))
end

_term_to_cube(term::Vector{<:SL.Atom}) = term
function _term_to_cube(term)
    term === SL.⊤ && return SL.Atom[]
    children = (term isa SL.Atom || term isa SL.Literal) ? (term,) : SL.grandchildren(term)
    return collect(SL.Atom, _lit_to_atom(c) for c in children)
end


# ---------------------------------------------------------------------------- #
#                    OFF-set construction (shared by seal + combine)           #
# ---------------------------------------------------------------------------- #
"""
    _collect_offset_cubes(buffers::Vector{_SeqBuffer}, exclude_idx::Int)
        -> Vector{Vector{SL.Atom}}

OFF-set esplicito per la classe `exclude_idx`: union dei cubi di TUTTE LE
ALTRE classi osservate finora in questo restart -- sia quelli già
sigillati in `phi_groups` (riconvertiti a cubo via `_term_to_cube`) sia
quelli ancora in `raw` (non ancora sigillati, ma già sicuramente
osservati e quindi confermati OFF per `exclude_idx`, essendo il modello
deterministico). Usata sia per sigillare un batch (`_seal_raw_batch!`) sia
per la ri-minimizzazione finale in fase di combine
(`_combine_phi_groups`).
"""
function _collect_offset_cubes(buffers::Vector{_SeqBuffer}, exclude_idx::Int)
    pieces = Vector{Vector{SL.Atom}}()
    for (j, other) in enumerate(buffers)
        j == exclude_idx && continue
        for grp in other.phi_groups
            isempty(grp) || append!(pieces, _term_to_cube.(grp))
        end
        isempty(other.raw) || append!(pieces, other.raw)
    end
    return pieces
end


# ---------------------------------------------------------------------------- #
#              STEP 1: sigillo isolato di un singolo batch di raw              #
# ---------------------------------------------------------------------------- #
"""
    _seal_raw_batch!(buf::_SeqBuffer, config::LumenConfig, M::Int;
                      offset::Vector{Vector{SL.Atom}}=Vector{Vector{SL.Atom}}()) -> Bool

STEP 1 della strategia seal-then-combine (vedi nota di patch in testa al
file). A differenza del vecchio `_fold_in!`, questa funzione NON combina
mai `buf.raw` con termini già sigillati in precedenza: minimizza
ESCLUSIVAMENTE `buf.raw` così com'è in questo momento (un singolo batch,
isolato), e il risultato -- una nuova φᵢ -- viene appeso in coda a
`buf.phi_groups`, senza toccare i gruppi già presenti.

    raw_batchᵢ → (Espresso, con offset dalle altre classi) → φᵢ
    push!(buf.phi_groups, φᵢ)
    empty!(buf.raw)

`offset` deve essere l'OFF-set esplicito costruito da
`_collect_offset_cubes` (o vuoto, se non servono garanzie contro le altre
classi); viene forwardato a `run_minimization(Val(:mitespresso), ...;
offset)`, che emette un PLA `.type fr`: `buf.raw` è l'ON-set confermato,
`offset` è l'OFF-set confermato, tutto il resto (incluse le righe di
questa STESSA classe non ancora enumerate, e i gruppi già sigillati di
questa stessa classe, che qui NON servono da offset perché sono ON, non
OFF) resta don't-care implicito.

Il valore di ritorno indica se la φᵢ appena sigillata sta, DA SOLA, sotto
il budget M -- NON è più (come nella vecchia versione) una garanzia sul
totale cumulativo della classe, perché con seal-then-combine il totale
cumulativo emerge solo dopo la combine finale (vedi nota di patch).
"""
function _seal_raw_batch!(
    buf::_SeqBuffer,
    config::LumenConfig,
    M::Int;
    offset::Vector{Vector{SL.Atom}}=Vector{Vector{SL.Atom}}()
)
    isempty(buf.raw) && return true

    n_raw = length(buf.raw)
    batch = copy(buf.raw)  # isolato: NON unito a nessun phi_groups precedente

    _dbg("seal_batch: raw_in_batch=", n_raw, " offset=", length(offset),
        " M=", M, " scheme=mitespresso, existing_phi_groups=", length(buf.phi_groups),
        " existing_sealed_terms=", _n_sealed_terms(buf))

    minimized = run_minimization(Val(:mitespresso), config, batch; offset=offset)

    minimized_cubes = convert(
        Vector{Vector{SL.Atom}},
        _term_to_cube.(collect(Any, minimized))
    )
    _verify_fold_coverage(batch, minimized_cubes; class_label=nothing, M=M)

    push!(buf.phi_groups, collect(Any, minimized))
    empty!(buf.raw)

    new_group_size = length(buf.phi_groups[end])
    _dbg("seal_batch: -> sealed new phi group, size=", new_group_size,
        " (", (new_group_size < M ? "OK, group under M" : "GROUP STILL OVER M"), ")",
        " total phi groups now=", length(buf.phi_groups),
        " cumulative sealed terms now=", _n_sealed_terms(buf))

    return new_group_size < M
end


# ---------------------------------------------------------------------------- #
#           STEP 2: combine finale dei φ sigillati + ri-minimizzazione         #
# ---------------------------------------------------------------------------- #
"""
    _combine_phi_groups(buffers::Vector{_SeqBuffer}, config::LumenConfig, M::Int;
                         reminimize::Bool=true) -> Vector{Vector{Any}}

STEP 2 della strategia seal-then-combine, eseguito UNA SOLA VOLTA per
restart, dopo che TUTTI i batch di TUTTE le classi sono stati sigillati
(sia quelli scattati durante l'enumerazione al raggiungimento di M, sia il
sigillo finale del raw residuo sotto M).

1. Per ogni classe, `combined = vcat(buf.phi_groups...)`: semplice unione
   letterale dei termini di tutte le φᵢ sigillate per quella classe.
   Nessuna minimizzazione avviene qui -- è pura concatenazione, quindi può
   solo essere ridondante (termini in più), mai lossy (non perde mai
   copertura, perché ogni φᵢ da sola copriva già il suo batch).

2. Se `reminimize=true` (default), viene fatta UNA sola passata aggiuntiva
   di `run_minimization` sul `combined` di ciascuna classe, con l'OFF-set
   esplicito costruito dai `combined` FINALI di tutte le altre classi (non
   più dai `phi_groups`/`raw` grezzi di pass, che a questo punto sono
   entrambi confluiti in `combined`). Questo squeeze finale rimuove la
   ridondanza fra batch diversi che il sigillo "a pezzi" ha
   necessariamente lasciato indietro (dato che ogni φᵢ non ha mai visto le
   altre).

   Se `reminimize=false`, viene restituito `combined` così com'è: resta
   comunque corretto (unione di minimizzazioni parziali, ciascuna sound),
   solo potenzialmente più grande del necessario -- utile se si vuole
   risparmiare l'ultima, potenzialmente costosa, chiamata a Espresso su un
   OFF-set grande quanto tutte le altre classi messe insieme.
"""
function _combine_phi_groups(
    buffers::Vector{_SeqBuffer},
    config::LumenConfig,
    M::Int;
    reminimize::Bool=true
)
    nclasses = length(buffers)
    combined_per_class = Vector{Vector{Any}}(undef, nclasses)
    @inbounds for ci in 1:nclasses
        combined_per_class[ci] = isempty(buffers[ci].phi_groups) ?
                                 Any[] : vcat(buffers[ci].phi_groups...)
    end

    _dbg("combine: pre_reminimize sizes=", length.(combined_per_class),
        " reminimize=", reminimize)

    reminimize || return combined_per_class

    final = Vector{Vector{Any}}(undef, nclasses)
    @inbounds for ci in 1:nclasses
        if isempty(combined_per_class[ci])
            final[ci] = combined_per_class[ci]
            continue
        end

        offset = Vector{Vector{SL.Atom}}()
        for (j, other) in enumerate(combined_per_class)
            j == ci && continue
            isempty(other) || append!(offset, _term_to_cube.(other))
        end

        cubes = convert(
            Vector{Vector{SL.Atom}},
            _term_to_cube.(combined_per_class[ci])
        )

        _dbg("combine: class_idx=", ci, " pre_reminimize_terms=", length(cubes),
            " offset=", length(offset))

        minimized = run_minimization(Val(:mitespresso), config, cubes; offset=offset)
        final[ci] = collect(Any, minimized)

        _dbg("combine: class_idx=", ci, " -> post_reminimize_terms=", length(final[ci]))
    end

    _dbg("combine: final sizes=", length.(final))
    return final
end


# ---------------------------------------------------------------------------- #
#         STEP 3: riconciliazione ESATTA contro il modello (nuovo)             #
# ---------------------------------------------------------------------------- #
"""
    _reconcile_against_model!(per_class_terms, config, model, ctx; max_apply_batch) -> Vector{Vector{Any}}

Perché serve: seal-then-combine (e allo stesso modo il vecchio fold
incrementale) può sigillare una φᵢ generalizzando su territorio non ancora
enumerato per NESSUNA classe in quel momento -- don't-care legittimo AL
MOMENTO del sigillo, ma che più avanti nello stesso restart può rivelarsi
appartenere davvero a un'altra classe. Quella φᵢ, una volta sigillata, non
viene mai ritoccata: il risultato osservabile è esattamente
`AVGSensitivity < 1` / `AVGSpecificity < 1` per alcune classi, perché le
regole finali non coincidono più esattamente con la funzione di decisione
del modello. Questo è INTRINSECO a qualunque strategia che sigilla/scarta
il raw prima che l'intero spazio sia stato enumerato -- non è specifico
del fold incrementale né di seal-then-combine.

Fix: dopo aver scelto il restart migliore, si fa UNA riconciliazione
esatta, in due sotto-passate sull'intero `ctx.n_total` (stesso ordine di
costo di UN restart, pagato una sola volta anziché N volte):

1. VERIFICA: per ogni cella dello spazio, calcola la classe vera (dal
   modello) e la confronta con quali regole finali la coprono. Ogni classe
   con anche una sola discrepanza viene marcata `dirty`.

2. RICOSTRUZIONE ESATTA (solo classi dirty): per ogni classe dirty si
   ricostruisce, in una seconda scansione, il suo ON-set ESATTO (tutte le
   celle la cui classe vera è quella) e si minimizza con OFF-set anch'esso
   ESATTO: per le altre classi dirty si usa il loro ON-set esatto appena
   ricostruito, per le classi NON dirty si usano le loro regole finali,
   che per definizione di "non dirty" sono già state verificate corrette
   su OGNI cella dello spazio al punto 1. Il risultato è un PLA `.type fr`
   senza don't-care residui per quella classe: Espresso è quindi
   matematicamente vincolato a produrre una copertura sound e completa,
   esattamente come farebbe Lumen minimizzando quella classe in un colpo
   solo.

   Nota deliberata: non si mescolano MAI i vecchi cubi-jolly (con
   wildcard) della classe dirty con minterm OFF specifici nella stessa
   chiamata -- sarebbe una specifica PLA contraddittoria (una riga jolly
   che "deve" coprire una cella dichiarata OFF). Per questo la classe
   dirty viene ricostruita da zero con dati esatti, non "corretta" a
   partire dalle sue vecchie φ.

Se nessuna classe risulta dirty, non viene fatta alcuna seconda passata né
alcuna chiamata extra a Espresso.
"""
function _reconcile_against_model!(
    per_class_terms::Vector{Vector{Any}},
    config::LumenConfig,
    model::SM.AbstractModel,
    ctx::NamedTuple;
    max_apply_batch::Int
)
    nclasses = length(ctx.classnames)
    class_index = Dict(c => i for (i, c) in enumerate(ctx.classnames))

    # ---- sotto-passata 1: verifica, individua le classi dirty ----
    dirty = falses(nclasses)

    i0 = 0
    while i0 < ctx.n_total
        this_chunk = clamp(min(max_apply_batch, ctx.n_total - i0), 1, ctx.n_total - i0)
        rows = Vector{NTuple{length(ctx.featurenames),eltype(ctx.thresholds[1])}}(undef, this_chunk)
        @inbounds for k in 1:this_chunk
            rows[k] = _combo_at(ctx.thrs_with_p, ctx.lens, ctx.strides, i0 + k - 1)
        end
        i0 += this_chunk

        tbl = NamedTuple{Tuple(ctx.featurenames)}(
            ntuple(j -> [r[j] for r in rows], length(ctx.featurenames))
        )
        d = PropositionalLogiset(tbl)
        preds = get_apply_function(config)(
            model, d;
            use_multithreads=get_use_multithreads(config),
            suppress_parity_warning=true
        )

        @inbounds for k in 1:this_chunk
            all(dirty) && break
            true_ci = get(class_index, preds[k], nothing)
            isnothing(true_ci) && continue

            truths_row = _truths_by_thresholds(rows[k], ctx.thresholds)
            cube = generate_disjunct(
                truths_row, ctx.thresholds, ctx.featurenames, ctx.op_families
            )

            for cj in 1:nclasses
                dirty[cj] && continue
                covered = !isempty(per_class_terms[cj]) && any(
                    t -> _cube_satisfies_term(cube, _term_to_cube(t)),
                    per_class_terms[cj]
                )
                if (cj == true_ci) != covered
                    dirty[cj] = true
                    _dbg("reconcile: class_idx=", cj, " marked DIRTY (cell mismatch found)")
                end
            end
        end
    end

    if !any(dirty)
        _dbg("reconcile: no dirty classes, model already exactly matched")
        return per_class_terms
    end

    _dbg("reconcile: dirty classes = ", ctx.classnames[dirty])

    # ---- sotto-passata 2: ON-set esatto per le sole classi dirty ----
    exact_on_dirty = Dict(ci => Vector{Vector{SL.Atom}}() for ci in findall(dirty))

    i0 = 0
    while i0 < ctx.n_total
        this_chunk = clamp(min(max_apply_batch, ctx.n_total - i0), 1, ctx.n_total - i0)
        rows = Vector{NTuple{length(ctx.featurenames),eltype(ctx.thresholds[1])}}(undef, this_chunk)
        @inbounds for k in 1:this_chunk
            rows[k] = _combo_at(ctx.thrs_with_p, ctx.lens, ctx.strides, i0 + k - 1)
        end
        i0 += this_chunk

        tbl = NamedTuple{Tuple(ctx.featurenames)}(
            ntuple(j -> [r[j] for r in rows], length(ctx.featurenames))
        )
        d = PropositionalLogiset(tbl)
        preds = get_apply_function(config)(
            model, d;
            use_multithreads=get_use_multithreads(config),
            suppress_parity_warning=true
        )

        @inbounds for k in 1:this_chunk
            true_ci = get(class_index, preds[k], nothing)
            (isnothing(true_ci) || !dirty[true_ci]) && continue

            truths_row = _truths_by_thresholds(rows[k], ctx.thresholds)
            cube = generate_disjunct(
                truths_row, ctx.thresholds, ctx.featurenames, ctx.op_families
            )
            push!(exact_on_dirty[true_ci], cube)
        end
    end

    # ---- ricostruzione: una minimizzazione esatta per ogni classe dirty ----
    for ci in findall(dirty)
        offset = Vector{Vector{SL.Atom}}()
        for cj in 1:nclasses
            cj == ci && continue
            if dirty[cj]
                append!(offset, exact_on_dirty[cj])
            elseif !isempty(per_class_terms[cj])
                append!(offset, _term_to_cube.(per_class_terms[cj]))
            end
        end

        _dbg("reconcile: rebuilding class_idx=", ci,
            " exact_on=", length(exact_on_dirty[ci]), " exact_offset=", length(offset))

        minimized = run_minimization(
            Val(:mitespresso), config, exact_on_dirty[ci]; offset=offset
        )
        per_class_terms[ci] = collect(Any, minimized)

        _dbg("reconcile: class_idx=", ci, " -> rebuilt exactly, size=",
            length(per_class_terms[ci]))
    end

    return per_class_terms
end


# ---------------------------------------------------------------------------- #
#                        one sequential pass (one restart)                     #
# ---------------------------------------------------------------------------- #
function _sequential_pass(
    config::LumenConfig,
    model::SM.AbstractModel,
    ctx::NamedTuple,
    seed::Integer;
    M::Int,
    max_apply_batch::Int,
    reminimize_after_combine::Bool=true
)
    nclasses = length(ctx.classnames)
    class_index = Dict(c => i for (i, c) in enumerate(ctx.classnames))
    buffers = [_SeqBuffer() for _ in 1:nclasses]

    perm = LazyPermutation(ctx.n_total; seed=seed)
    ok = true

    _dbg("pass start: seed=", seed, " n_total=", ctx.n_total, " M=", M,
        " max_apply_batch=", max_apply_batch,
        " reminimize_after_combine=", reminimize_after_combine)

    n_seals = 0
    i0 = 0
    while i0 < ctx.n_total
        # Il "room" per il prossimo chunk dipende solo dal raw non ancora
        # sigillato: i phi_groups già sigillati non contano più verso M,
        # dato che non vengono più ri-processati insieme a nuovo raw.
        room = M - maximum(_raw_occupancy(b) for b in buffers; init=0)
        this_chunk = clamp(min(room, max_apply_batch, ctx.n_total - i0), 1, ctx.n_total - i0)

        rows = Vector{NTuple{length(ctx.featurenames),eltype(ctx.thresholds[1])}}(undef, this_chunk)
        @inbounds for k in 1:this_chunk
            linidx = permute(perm, i0 + k - 1)
            rows[k] = _combo_at(ctx.thrs_with_p, ctx.lens, ctx.strides, linidx)
        end
        i0 += this_chunk

        tbl = NamedTuple{Tuple(ctx.featurenames)}(
            ntuple(j -> [r[j] for r in rows], length(ctx.featurenames))
        )
        d = PropositionalLogiset(tbl)
        preds = get_apply_function(config)(
            model, d;
            use_multithreads=get_use_multithreads(config),
            suppress_parity_warning=true
        )

        @inbounds for k in 1:this_chunk
            ci = get(class_index, preds[k], nothing)
            isnothing(ci) && continue

            truths_row = _truths_by_thresholds(rows[k], ctx.thresholds)
            cube = generate_disjunct(
                truths_row, ctx.thresholds, ctx.featurenames, ctx.op_families
            )
            buf = buffers[ci]
            push!(buf.raw, cube)
            buf.n_rows_seen += 1

            if _raw_occupancy(buf) >= M
                _dbg("trigger seal: class=", ctx.classnames[ci],
                    " raw_occupancy=", _raw_occupancy(buf), " (i0=", i0, "/", ctx.n_total, ")")
                n_seals += 1
                # OFF-set esplicito dalle altre classi, costruito PRIMA di
                # sigillare questo batch, così Espresso non deve mai
                # indovinare per complemento assoluto.
                offset = _collect_offset_cubes(buffers, ci)
                ok &= _seal_raw_batch!(buf, config, M; offset=offset)
                if !ok
                    _dbg("restart seed=", seed, " ABORTED: class=", ctx.classnames[ci],
                        " single sealed batch did not fit under M=", M)
                    combined_partial = _combine_phi_groups(
                        buffers, config, M; reminimize=false
                    )
                    return combined_partial, false
                end
            end
        end
    end

    # Sigillo finale: raw residuo sotto M per ogni classe, che non aveva
    # ancora fatto scattare il trigger `>= M` durante l'enumerazione.
    for (ci, buf) in enumerate(buffers)
        offset = _collect_offset_cubes(buffers, ci)
        ok &= _seal_raw_batch!(buf, config, M; offset=offset)
        if !ok
            _dbg("restart seed=", seed, " ABORTED at final seal: class=",
                ctx.classnames[ci], " did not fit under M")
        end
    end

    # STEP 2 (unico, a fine pass): combina tutti i φ sigillati di ogni
    # classe e, se richiesto, fa la ri-minimizzazione finale.
    per_class_terms = _combine_phi_groups(
        buffers, config, M; reminimize=reminimize_after_combine
    )

    _dbg("pass end: seed=", seed, " ok=", ok, " n_seals=", n_seals,
        " per_class_sizes=", length.(per_class_terms),
        " n_rows_seen=", [b.n_rows_seen for b in buffers])

    return per_class_terms, ok
end


# ---------------------------------------------------------------------------- #
#                                public entry point                            #
# ---------------------------------------------------------------------------- #
"""
    super_lumen(config::LumenConfig, model::SM.AbstractModel;
                M::Int=100_000,
                N::Int=10,
                max_apply_batch::Int=min(M, 4096),
                score::Symbol=:n_terms,
                reminimize_after_combine::Bool=true,
                reconcile_final::Bool=true) -> SM.DecisionSet

Vedi la nota di patch a inizio file per la differenza rispetto alla
revisione precedente: ogni batch di raw cube viene ora minimizzato in
isolamento e "sigillato" a parte (`raw_batchᵢ → φᵢ`, mai
`φ + raw_batchᵢ → φ`), e solo a fine restart tutte le φᵢ di ciascuna
classe vengono unite e, se `reminimize_after_combine=true`, ri-minimizzate
una volta sola.

`reminimize_after_combine`: se `false`, la formula finale di ogni classe è
semplicemente l'unione letterale (`vcat`) di tutte le φ sigillate durante
il restart -- corretta ma potenzialmente più grande del necessario, e più
economica (nessuna chiamata extra a Espresso a fine pass). Se `true`
(default), viene eseguita un'ultima chiamata a Espresso per classe per
comprimere quell'unione.

`reconcile_final`: se `true` (default), dopo aver scelto il restart
migliore viene eseguita `_reconcile_against_model!`, che verifica le
regole finali contro il modello reale su TUTTO lo spazio enumerato e
ricostruisce esattamente (senza approssimazione) le classi che
presentassero anche una sola discrepanza -- questo è ciò che garantisce
`AVGSensitivity == AVGSpecificity == 1.0`, allineando SuperLumen a Lumen
sulla correttezza. Costa una scansione completa di `ctx.n_total` (più una
seconda, ristretta alle sole classi "dirty"), pagata UNA VOLTA a fine
pipeline, non per ogni restart. Disattivarlo (`false`) salta questo
controllo e riporta il comportamento "best-effort" della revisione
precedente (più veloce, ma senza garanzia di sensitivity/specificity=1).

Tutto il resto -- restart loop, scoring, DecisionSet assembly -- è
invariato rispetto alla revisione precedente.
"""
function super_lumen(
    config::LumenConfig,
    model::SM.AbstractModel;
    M::Int=100_000,
    N::Int=10,
    max_apply_batch::Int=min(M, 4096),
    score::Symbol=:n_terms,
    reminimize_after_combine::Bool=true,
    reconcile_final::Bool=true
)
    score === :n_terms || throw(ArgumentError(
        "Only score=:n_terms is currently supported."
    ))
    M >= 1 || throw(ArgumentError("M must be positive, got $M"))
    N >= 1 || throw(ArgumentError("N must be positive, got $N"))

    ctx = _prepare_sequential_context(config, model)

    _dbg("super_lumen start: M=", M, " N=", N,
        " max_apply_batch=", max_apply_batch, " scheme=mitespresso (hardcoded)",
        " reminimize_after_combine=", reminimize_after_combine,
        " reconcile_final=", reconcile_final)

    best_terms = nothing
    best_score = typemax(Int)
    n_discarded = 0

    for r in 1:N
        seed = rand(UInt64)
        per_class_terms, ok = _sequential_pass(
            config, model, ctx, seed;
            M, max_apply_batch, reminimize_after_combine
        )
        if !ok
            n_discarded += 1
            _dbg("restart ", r, "/", N, " (seed=", seed, ") DISCARDED (over budget)")
            continue
        end

        this_score = sum(length, per_class_terms)
        _dbg("restart ", r, "/", N, " (seed=", seed, ") score=", this_score,
            " best_so_far=", best_score)
        if this_score < best_score
            best_score = this_score
            best_terms = per_class_terms
            _dbg("  -> new best (score=", this_score, ")")
        end
    end

    _dbg("super_lumen: ", N - n_discarded, "/", N, " restarts kept, best_score=", best_score)

    isnothing(best_terms) && error(
        "super_lumen: all $N restarts exceeded the memory budget " *
        "(M=$M) even after sealing each batch in isolation. Raise M, " *
        "or increase N to give more random orderings a chance to seal " *
        "each batch under budget."
    )

    if reconcile_final
        _dbg("super_lumen: running final reconciliation against the model " *
             "(guarantees sensitivity=specificity=1 for reconstructed classes)")
        best_terms = _reconcile_against_model!(
            best_terms, config, model, ctx; max_apply_batch
        )
        best_score = sum(length, best_terms)
        _dbg("super_lumen: post-reconciliation score=", best_score)
    end

    valid_mask = .!isempty.(best_terms)
    classes = ctx.classnames[valid_mask]

    float_type = get_float_type(config)
    formulas = Vector{Vector{Union{
        SL.LeftmostConjunctiveForm{SL.Atom{float_type}},
        SyntaxStructure
    }}}(best_terms[valid_mask])

    _dbg("super_lumen done: nclasses_kept=", length(classes),
        " total_terms=", best_score)

    return SM.DecisionSet(
        SM.Rule.(SL.LeftmostDisjunctiveForm.(formulas), classes)
    )
end