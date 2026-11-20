# Run explicitly once. No package installation is hidden inside an analysis.
needed<-c('lpSolve','quadprog','factorstochvol','posterior')
missing<-needed[!vapply(needed,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing))install.packages(missing,repos='https://cloud.r-project.org')
stopifnot(all(vapply(needed,requireNamespace,logical(1),quietly=TRUE)))
print(data.frame(package=needed,version=vapply(needed,function(x)as.character(packageVersion(x)),character(1))))
