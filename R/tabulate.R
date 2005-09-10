cl_tabulate <-
function(x)
{
    values <- unique(x)
    counts <- tabulate(match(x, values))
    data.frame(values = I(values), counts = counts)
}
