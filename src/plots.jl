function plotresults(parametervec, xlabel; truevalue = [])

  if length(truevalue) > 0
    p = plot(x = parametervec, Geom.density,
     xintercept = truevalue,
     Geom.vline(color = "red", size=1mm),
     Guide.xlabel(xlabel))
   else
     p = plot(x = parametervec, Geom.density,
      Guide.xlabel(xlabel))
   end

  return p
end


function plotmodelposterior(res::ABCSMCmodelresults)
    p = plot(res.ModelProb, x=:Model, y = :Probability, Geom.bar,
    Theme(bar_spacing = 0.2cm,
    default_color = RGBA(0.5, 0.5, 0.5, 0.8),
    major_label_font_size = 16pt,
    minor_label_font_size = 12pt))

    return p
end

function plotparameterposterior(res::ABCSMCmodelresults, model = 1)

    DF = stack(res.Posterior[model].Parameters)
    plot(DF, x=:value, xgroup=:variable,
    Geom.subplot_grid(Geom.histogram(bincount = 30), free_x_axis=true),
    Theme(
    default_color = RGBA(0.5, 0.5, 0.5, 0.9),
    major_label_font_size = 12pt,
    minor_label_font_size = 8pt))

end


function plotparameterposterior(res::ABCSMCresults)

    DF = stack(DataFrame(res.parameters))
    plot(DF, x=:value, xgroup=:variable,
    Geom.subplot_grid(Geom.histogram(bincount = 30), free_x_axis=true),
    Theme(
    default_color = RGBA(0.5, 0.5, 0.5, 0.9),
    major_label_font_size = 12pt,
    minor_label_font_size = 8pt))

end
