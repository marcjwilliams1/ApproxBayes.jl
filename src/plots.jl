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
