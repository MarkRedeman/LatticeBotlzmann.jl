using Documenter, lbm

makedocs(
    modules = [lbm],
    format = :html,
    checkdocs = :exports,
    sitename = "lbm.jl",
    pages = Any["index.md"],
)

deploydocs(repo = "github.com/MarkRedeman/lbm.jl.git")
