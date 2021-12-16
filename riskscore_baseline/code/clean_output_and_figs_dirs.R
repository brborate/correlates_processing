file.remove(paste0("output/", grep(list.files(path="output"), pattern='.RData', invert=TRUE, value=TRUE)))
file.remove((list.files(here("figs"), full.names = TRUE)))


