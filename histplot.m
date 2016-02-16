function histplotout = histplot(values,range)

histc(values,range)
histplotout = histc(values,range)/numel(values);


end

