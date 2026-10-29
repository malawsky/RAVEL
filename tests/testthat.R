library(testthat)
library(RAVEL)

# Match the existing development runner before tests call lavaan directly.
# Restricted hosts can report NA from detectCores(), poisoning lavaan's cached
# default options if its first initialization precedes ravel_fit().
RAVEL:::.ensure_lavaan_runtime()

test_check("RAVEL")
