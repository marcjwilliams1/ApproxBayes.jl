"""
    plotmodelposterior(results; <keyword arguments>)

Plot the posterior probabalities of each model.
...
## Arguments
- `save = false`: Set to `true` if you want the plot to be saved
- `dir = ""`: Directory where the plot will be saved to. Default is the current working directory.
- `plotname = ""`: Name to call plot when saving, default depends on the algorithm used.
...
"""
function plotmodelposterior(res::ABCSMCmodelresults; save = false, dir = "", plotname = "ABCSMCmodelposteriors")

    #pyplot()

    DF = DataFrame(Model = map(x -> "$x", 1:length(res.modelprob)), Probability = res.modelprob)

    Plots.bar(DF[:Model], DF[:Probability],
    title="Model Probabilities", yaxis = ("Probability"),
    linecolor = :white, fillcolor = RGBA(0.5, 0.5, 0.5, 0.8),
    markerstrokecolor=:white, titlefont = font(12, "Calibri"), ytickfont = font(10, "Calibri"), xtickfont = font(10, "Calibri"), legend = false, grid = false)

    if save == true
        Plots.savefig(joinpath(dir, "$(plotname).pdf"))
    end

    return p
end

function plotmodelposterior(res::ABCrejectionmodelresults; save = false, dir = "", plotname = "ABCRejectionmodelposteriors")

    #pyplot()
    DF = DataFrame(Model = map(x -> "$x", 1:length(res.modelfreq)), Probability = res.modelfreq)

    Plots.bar(DF[:Model], DF[:Probability],
    title="Model Probabilities", yaxis = ("Probability"),
    linecolor = :white, fillcolor = RGBA(0.5, 0.5, 0.5, 0.8),
    markerstrokecolor=:white, titlefont = font(12, "Calibri"), ytickfont = font(10, "Calibri"), xtickfont = font(10, "Calibri"), legend = false, grid = false)

    if save == true
        Plots.savefig(joinpath(dir, "$(plotname).pdf"))
    end

    return p
end

"""
    plotparameterposterior(results [,model = 1]; <keyword arguments>)

Plot the parameter posterior values. If algorithm is a model selection algorithm, model needs to be specified.
...
## Arguments
- `save = false`: Set to `true` if you want the plot to be saved
- `dir = ""`: Directory where the plot will be saved to. Default is the current working directory.
- `plotname = ""`: Name to call plot when saving, default depends on the algorithm used.
...
"""
function plotparameterposterior(res::ABCrejectionmodelresults, model = 1; save = false, dir = "", plotname = "ABCRejectionparameterposteriors")
    #pyplot()
    nparams = size(res.parameters[model])[2]

    Plots.histogram(Array(res.parameters[model][:, 1:nparams]), nbins = 20, layout = nparams, normed = true,
    title=hcat(map(x -> "Parameter $x", 1:nparams)...),
    linecolor = :white, fillcolor = RGBA(0.75, 0.3, 0.3),
    markerstrokecolor=:white, titlefont = font(12, "Calibri"), ytickfont = font(10, "Calibri"), xtickfont = font(10, "Calibri"), legend = false)

    if save == true
        Plots.savefig(joinpath(dir, "$(plotname)-model$(model).pdf"))
    end
end

function plotparameterposterior(res::ABCrejectionresults; save = false, dir = "", plotname = "ABCRejectionparameterposteriors")
    nparams = size(res.parameters)[2]
    #pyplot()
    Plots.histogram(Array(res.parameters[:, 1:nparams]),
    nbins = 20, layout = nparams, normed = true,
    title=hcat(map(x -> "Parameter $x", 1:nparams)...),
    linecolor = :white, fillcolor = RGBA(0.75, 0.3, 0.3),
    markerstrokecolor=:white, titlefont = font(12, "Calibri"), ytickfont = font(10, "Calibri"), xtickfont = font(10, "Calibri"), legend = false)

    if save == true
        Plots.savefig(joinpath(dir, "$(plotname).pdf"))
    end
end

function plotparameterposterior(res::ABCSMCmodelresults, model = 1; save = false, dir = "", plotname = "ABCSMCparameterposteriors")
    #pyplot()
    nparams = size(res.parameters[model])[2]

    Plots.histogram(Array(res.parameters[model])[:, 1:nparams], nbins = 20, layout = nparams, weights = res.weights[model],
    title=hcat(map(x -> "Parameter $x", 1:nparams)...),
    linecolor = :white, fillcolor = RGBA(0.75, 0.3, 0.3),
    markerstrokecolor=:white, titlefont = font(12, "Calibri"), ytickfont = font(10, "Calibri"), xtickfont = font(10, "Calibri"), legend = false)

    if save == true
        Plots.savefig(joinpath(dir, "$(plotname)-model$(model).pdf"))
    end
end

function plotparameterposterior(res::ABCSMCresults; save = false, dir = "", plotname = "ABCSMCparameterposteriors")
    #pyplot()
    nparams = size(res.parameters)[2]

    Plots.histogram(Array(res.parameters)[:, 1:nparams], nbins = 20, layout = nparams, weights = res.weights,
    title=hcat(map(x -> "Parameter $x", 1:nparams)...),
    linecolor = :white, fillcolor = RGBA(0.75, 0.3, 0.3),
    markerstrokecolor=:white, titlefont = font(12, "Calibri"), ytickfont = font(10, "Calibri"), xtickfont = font(10, "Calibri"), legend = false)

    if save == true
        Plots.savefig(joinpath(dir, "$(plotname).pdf"))
    end
end
