function plotresults(parametervec, xlabel, truevalue)

  p = plot(x = parametervec, Geom.density,
   xintercept = [truevalue],
   Geom.vline(color = "red", size=1mm),
   Guide.xlabel(xlabel))

  return p
end
