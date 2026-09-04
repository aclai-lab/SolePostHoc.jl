using SoleData
const SD = SoleData

using SoleData.Artifacts

# fill your Artifacts.toml file;
fillartifacts() # comment this line in debug

espressobinary = try
    joinpath(SD.load(SD.MITESPRESSOLoader()), "espresso")
catch e
    error("Failed to setup espresso binary: $e")
end