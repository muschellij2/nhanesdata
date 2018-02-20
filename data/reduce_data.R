
objects = sub(".rda",  "",
        list.files("data", pattern = ".rda"))
iobj = objects[length(objects)]

for (iobj in objects) {
  fname = file.path("data", paste0(iobj, ".rda"))
  e = environment()
  load(fname, envir = e)
  df = e[[iobj]]
  for (icol in colnames(df)) {
    x = df[, icol]
    if (is.numeric(x)) {
      int_x = as.integer(x)
      if (all(x == int_x, na.rm = TRUE)) {
        df[, icol] = int_x
      } else {
        print(icol)
      }
    }
  }
  save(list = iobj, envir = e, file = fname, compress = "xz")
  # save()
}

