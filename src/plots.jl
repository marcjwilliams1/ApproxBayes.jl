function plotresults(parametervec, xlabel; truevalue = [])

  if length(truevalue) > 0
    p = Gadfly.plot(x = parametervec, Geom.density,
     xintercept = truevalue,
     Geom.vline(color = "red", size=1mm),
     Guide.xlabel(xlabel))
   else
     p = Gadfly.plot(x = parametervec, Geom.density,
      Guide.xlabel(xlabel))
   end

  return p
end


function plotmodelposterior(res::ABCSMCmodelresults)
    p = Gadfly.plot(res.ModelProb, x=:Model, y = :Probability, Geom.bar,
    Theme(bar_spacing = 0.2cm,
    default_color = RGBA(0.5, 0.5, 0.5, 0.8),
    major_label_font_size = 16pt,
    minor_label_font_size = 12pt))

    return p
end

function plotparameterposterior(res::ABCSMCmodelresults, model = 1)

    nparams = size(res.Posterior[model].Parameters)[2]

    Plots.histogram(Array(res.Posterior[model].Parameters)[:, 1:nparams], nbins = 20, layout = nparams,
    title=hcat(map(x -> "Parameter $x", 1:nparams)...),
    linecolor = :white, fillcolor = RGBA(0.75, 0.3, 0.3),
    markerstrokecolor=:white, titlefont = font(12, "Calibri"), ytickfont = font(10, "Calibri"), xtickfont = font(10, "Calibri"), legend = false)

end

function plotparameterposterior(res::ABCSMCresults)

    nparams = size(res.parameters)[2]

    Plots.histogram(Array(res.parameters)[:, 1:nparams], nbins = 20, layout = nparams,
    title=hcat(map(x -> "Parameter $x", 1:nparams)...),
    linecolor = :white, fillcolor = RGBA(0.75, 0.3, 0.3),
    markerstrokecolor=:white, titlefont = font(12, "Calibri"), ytickfont = font(10, "Calibri"), xtickfont = font(10, "Calibri"), legend = false)

end
