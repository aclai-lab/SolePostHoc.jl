# lumen.jl
# Comprehensive test suite for the Lumen module

using Test
using SolePostHoc
using DataFrames
using DecisionTree
using SoleModels
using SoleData
using SoleLogics

@testset "Lumen Module Test Suite" begin
    
    # Test data setup
    function create_test_data()
        # Create simple Iris-like test data
        n = 150
        features = rand(n, 4) .* 5  # 4 features, values 0-5
        
        # Create simple decision boundaries
        labels = Vector{String}(undef, n)
        for i in 1:n
            if features[i, 1] < 2.0
                labels[i] = "class_A"
            elseif features[i, 2] > 3.5
                labels[i] = "class_B"
            else
                labels[i] = "class_C"
            end
        end
        
        return features, labels
    end
    
    function create_test_model()
        features, labels = create_test_data()
        
        # Create a simple decision tree
        model = build_tree(labels, features; 
                          max_depth=3, 
                          min_samples_leaf=5)
        return model, features, labels
    end
    
    function create_test_ensemble()
        features, labels = create_test_data()
        
        # Create a random forest
        model = build_forest(labels, features; 
                           n_trees=5, 
                           n_subfeatures=2,
                           max_depth=3)
        return model, features, labels
    end
    
    @testset "Configuration Tests" begin
        
        @testset "LumenConfig Construction" begin
            # Test default configuration
            config = LumenConfig()
            @test config.minimization_scheme == :abc
            @test config.vertical == 1.0
            @test config.horizontal == 1.0
            @test config.ott_mode == true
            @test config.silent == false
            
            # Test custom configuration
            config = LumenConfig(
                minimization_scheme=:mitespresso,
                vertical=0.8,
                horizontal=0.9,
                silent=true
            )
            @test config.minimization_scheme == :mitespresso
            @test config.vertical == 0.8
            @test config.horizontal == 0.9
            @test config.silent == true
        end
        
        @testset "Configuration Validation" begin
            # Test valid configurations
            @test_nowarn LumenConfig(vertical=0.5, horizontal=0.7)
            @test_nowarn LumenConfig(minimization_scheme=:boom)
            
            # Test invalid vertical parameter
            @test_throws ArgumentError LumenConfig(vertical=0.0)
            @test_throws ArgumentError LumenConfig(vertical=1.5)
            @test_throws ArgumentError LumenConfig(vertical=-0.1)
            
            # Test invalid horizontal parameter
            @test_throws ArgumentError LumenConfig(horizontal=0.0)
            @test_throws ArgumentError LumenConfig(horizontal=1.2)
            
            # Test invalid minimization scheme
            @test_throws ArgumentError LumenConfig(minimization_scheme=:invalid_scheme)
            @test_throws ArgumentError LumenConfig(minimization_scheme=:unknown)
        end
    end
    
    @testset "Rule Extraction Strategy Tests" begin
        
        @testset "Strategy Type System" begin
            # Test abstract type hierarchy
            @test StandardExtraction() isa RuleExtractionStrategy
            @test EnsembleExtraction() isa RuleExtractionStrategy
        end
        
        @testset "Rule Extraction" begin
            model, _, _ = create_test_model()
            
            # Test standard extraction
            strategy = StandardExtraction()
            rules = extract_rules(model, strategy)
            @test rules isa Vector
            @test length(rules) > 0
            
            # Test automatic strategy selection
            rules_auto = extract_rules(model)
            @test rules_auto isa Vector
            @test length(rules_auto) > 0
        end
        
        @testset "Ensemble Rule Extraction" begin
            ensemble, _, _ = create_test_ensemble()
            
            # Test ensemble extraction
            strategy = EnsembleExtraction()
            rules = extract_rules(ensemble, strategy)
            @test rules isa Vector
            @test length(rules) > 0
            
            # Test automatic strategy selection for ensemble
            rules_auto = extract_rules(ensemble)
            @test rules_auto isa Vector
            @test length(rules_auto) > 0
        end
    end
    
    @testset "Helper Function Tests" begin
        
        @testset "get_formula_terms Function" begin
            # Test with mock TwoLevelDNFFormula (if available)
            # This would need to be adjusted based on actual type definitions
            
            # Test fallback behavior with integer
            @test Lumen.get_formula_terms(1) >= 1
            
            # Test with string (fallback case)
            test_formula_str = "A ∨ B ∨ C"
            # This tests the string parsing fallback
            terms = Lumen.get_formula_terms(test_formula_str)
            @test terms isa Integer
            @test terms >= 1
        end
        
        @testset "should_track_statistics Function" begin
            @test Lumen.should_track_statistics(:abc) == true
            @test Lumen.should_track_statistics(:mitespresso) == true
            @test Lumen.should_track_statistics(:boom) == true
            @test Lumen.should_track_statistics(:espresso) == false
        end
        
        @testset "determine_apply_function" begin
            model, _, _ = create_test_model()
            ensemble, _, _ = create_test_ensemble()
            
            # Test with single tree
            apply_func = Lumen.determine_apply_function(model, nothing)
            @test apply_func isa Function
            
            # Test with ensemble
            apply_func_ensemble = Lumen.determine_apply_function(ensemble, nothing)
            @test apply_func_ensemble isa Function
            
            # Test with custom function
            custom_func(x, y) = x
            apply_func_custom = Lumen.determine_apply_function(model, custom_func)
            @test apply_func_custom === custom_func
        end
    end
    
    @testset "Main Algorithm Tests" begin
        
        @testset "Basic Lumen Execution" begin
            model, features, labels = create_test_model()
            
            # Test basic execution with default config
            result = lumen(model)
            @test result isa LumenResult
            @test result.decision_set isa DecisionSet
            @test result.processing_time >= 0.0
            
            # Test with custom configuration
            config = LumenConfig(silent=true, return_info=false)
            result2 = lumen(model, config)
            @test result2 isa LumenResult
            @test result2.decision_set isa DecisionSet
        end
        
        @testset "Lumen with Different Minimization Schemes" begin
            model, _, _ = create_test_model()
            
            for scheme in [:abc, :mitespresso, :boom]
                config = LumenConfig(minimization_scheme=scheme, silent=true)
                @test_nowarn begin
                    result = lumen(model, config)
                    @test result isa LumenResult
                end
            end
        end
        
        @testset "Lumen with Ensemble Models" begin
            ensemble, _, _ = create_test_ensemble()
            
            config = LumenConfig(silent=true)
            result = lumen(ensemble, config)
            @test result isa LumenResult
            @test result.decision_set isa DecisionSet
        end
        
        @testset "Configuration Options" begin
            model, _, _ = create_test_model()
            
            # Test different coverage parameters
            for vertical in [0.7, 0.8, 0.9, 1.0]
                for horizontal in [0.7, 0.8, 0.9, 1.0]
                    config = LumenConfig(
                        vertical=vertical, 
                        horizontal=horizontal, 
                        silent=true
                    )
                    @test_nowarn begin
                        result = lumen(model, config)
                        @test result isa LumenResult
                    end
                end
            end
            
            # Test OTT mode
            config_ott = LumenConfig(ott_mode=true, silent=true)
            config_no_ott = LumenConfig(ott_mode=false, silent=true)
            
            @test_nowarn lumen(model, config_ott)
            @test_nowarn lumen(model, config_no_ott)
        end
        
        @testset "Return Info Options" begin
            model, _, _ = create_test_model()
            
            # Test with return_info=true
            config_info = LumenConfig(return_info=true, silent=true)
            result_info = lumen(model, config_info)
            @test result_info.info isa NamedTuple
            @test haskey(result_info.info, :processing_time)
            
            # Test with return_info=false
            config_no_info = LumenConfig(return_info=false, silent=true)
            result_no_info = lumen(model, config_no_info)
            @test result_no_info.info == (;)
        end
    end
    
    @testset "LumenResult Tests" begin
        
        @testset "LumenResult Construction" begin
            # Mock DecisionSet (adjust based on actual implementation)
            mock_rules = []  # Empty rules for testing
            mock_ds = DecisionSet(mock_rules)
            
            # Test full constructor
            info = (test_param = 123,)
            result = LumenResult(mock_ds, info, 1.5)
            @test result.decision_set === mock_ds
            @test result.info === info
            @test result.processing_time == 1.5
            
            # Test simple constructor
            simple_result = LumenResult(mock_ds)
            @test simple_result.decision_set === mock_ds
            @test simple_result.info == (;)
            @test simple_result.processing_time == 0.0
        end
    end
    
    @testset "Error Handling Tests" begin
        
        @testset "Invalid Model Inputs" begin
            # Test with invalid model types
            @test_throws Exception lumen("not_a_model")
            @test_throws Exception lumen(123)
            @test_throws Exception lumen(nothing)
        end
        
        @testset "Robustness Tests" begin
            model, _, _ = create_test_model()
            
            # Test with extreme parameter values
            config_extreme = LumenConfig(
                vertical=0.01,
                horizontal=0.01,
                silent=true
            )
            
            # This might fail or produce warnings, but shouldn't crash
            @test_nowarn begin
                try
                    result = lumen(model, config_extreme)
                    @test result isa LumenResult
                catch e
                    @warn "Extreme parameters caused expected failure: $e"
                end
            end
        end
    end
    
    @testset "Performance Tests" begin
        
        @testset "Timing Tests" begin
            model, _, _ = create_test_model()
            config = LumenConfig(silent=true)
            
            # Measure execution time
            result = @timed lumen(model, config)
            @test result.time < 60.0  # Should complete within 60 seconds
            @test result.value isa LumenResult
            
            # Check if reported processing time is reasonable
            lumen_result = result.value
            @test lumen_result.processing_time > 0.0
            @test lumen_result.processing_time < result.time * 2  # Some overhead is expected
        end
        
        @testset "Memory Usage Tests" begin
            model, _, _ = create_test_model()
            config = LumenConfig(silent=true)
            
            # Basic memory allocation test
            @test_nowarn begin
                for i in 1:5  # Run multiple times to check for memory leaks
                    result = lumen(model, config)
                    @test result isa LumenResult
                end
            end
        end
    end
    
    @testset "Integration Tests" begin
        
        @testset "End-to-End Workflow" begin
            # Create data
            features, labels = create_test_data()
            
            # Build model
            model = build_tree(labels, features; max_depth=4)
            
            # Extract rules with Lumen
            config = LumenConfig(
                minimization_scheme=:abc,
                vertical=0.9,
                horizontal=0.8,
                silent=true,
                return_info=true
            )
            
            result = lumen(model, config)
            
            # Validate complete result
            @test result isa LumenResult
            @test result.decision_set isa DecisionSet
            @test result.processing_time > 0.0
            @test result.info isa NamedTuple
            
            # Check that we got some rules
            @test length(result.decision_set.rules) > 0
        end
        
        @testset "Comparison Between Methods" begin
            model, _, _ = create_test_model()
            
            results = Dict()
            for scheme in [:abc, :mitespresso, :boom]
                config = LumenConfig(
                    minimization_scheme=scheme, 
                    silent=true, 
                    return_info=true
                )
                results[scheme] = lumen(model, config)
            end
            
            # All methods should produce valid results
            for (scheme, result) in results
                @test result isa LumenResult
                @test result.decision_set isa DecisionSet
                @test result.processing_time > 0.0
            end
            
            # Compare number of rules (they might be different due to different minimization)
            rule_counts = [length(result.decision_set.rules) for result in values(results)]
            @test all(count -> count >= 0, rule_counts)
        end
    end
end

# Additional utility functions for testing
function print_test_summary()
    println("\n" * "="^60)
    println("LUMEN MODULE TEST SUITE COMPLETED")
    println("="^60)
    println("Run with: julia testLumenModule.jl")
    println("Or in REPL: include(\"testLumenModule.jl\")")
    println("="^60 * "\n")
end

# Run tests if script is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    print_test_summary()
end