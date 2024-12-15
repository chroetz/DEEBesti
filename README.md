# DEEBesti

R implementations of methods for learning differential equations. This package is typically used via the command line interface provide by the package DEEBcmd.

Alternatively, a typicall useage of this package looks like this:

```r
DEEBesti::runOne(
    dbPath = "path/to/deebDb", # Path a DEEB database folder
    truthNrFilter = 1:10, # process only these truth numbers
    obsNr = 2:3,  # process only these observation numbers
    model = "lorenz63", # process this statistical model / dynamical system
    method = "Node1", # process this method Opts file from the DEEB database
    expansionNr = 1 # if the method file contains multiple methods (e.g., different hyperparameter settings), this identifies which method to run
))
```

