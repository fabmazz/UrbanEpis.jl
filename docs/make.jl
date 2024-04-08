using UrbanEpis
using Documenter

DocMeta.setdocmeta!(UrbanEpis, :DocTestSetup, :(using UrbanEpis); recursive=true)

makedocs(;
    modules=[UrbanEpis],
    authors="Fabio Mazza",
    sitename="UrbanEpis.jl",
    format=Documenter.HTML(;
        canonical="https://fabmazza.github.io/UrbanEpis.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/fabmazza/UrbanEpis.jl",
    devbranch="main",
)
