using DecayAngles
using Documenter

DocMeta.setdocmeta!(DecayAngles, :DocTestSetup, :(using DecayAngles); recursive = true)

const page_rename = Dict("developer.md" => "Developer docs") # Without the numbers
const numbered_pages = [
    file for file in readdir(joinpath(@__DIR__, "src")) if
    file != "index.md" && splitext(file)[2] == ".md"
]

makedocs(;
    modules = [DecayAngles],
    authors = "Mikhail Mikhasenko <mikhail.mikhasenko@cern.ch>",
    repo = "https://github.com/mmikhasenko/DecayAngles.jl/blob/{commit}{path}#{line}",
    sitename = "DecayAngles.jl",
    format = Documenter.HTML(; canonical = "https://mmikhasenko.github.io/DecayAngles.jl"),
    pages = ["index.md"; numbered_pages],
)

deploydocs(; repo = "github.com/mmikhasenko/DecayAngles.jl")
