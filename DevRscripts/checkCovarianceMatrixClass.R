library(ribModel)

cm <- new (CovarianceMatrix)
A = matrix(c(25,15,-5, 15,18,0, -5,0,11), nrow = 3)
cm$setCovarianceMatrix(A)
cm$printCovarianceMatrix()
cm$choleskiDecomposition()
cm$printCholeskiMatrix()