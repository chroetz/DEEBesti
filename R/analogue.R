predictDirectDeriv <- function(states, parms, hyperParms) {
  t(apply(states, 1, \(state) {
    analogue <-
      extendAnalogue(
        makeTrajs(0, matrix(state, nrow=1)),
        parms,
        requireSteps = hyperParms$derivOrder)
    analogue <- analogue[seq_len(hyperParms$derivOrder+1), ]
    polyInterpCoeffs <- polynomialInterpolation(analogue$time, analogue$state)
    polyInterpCoeffs[2,] # derivative at 0 of polynomial is linear coefficient (second coeff)
  }))
}


extendAnalogue <- function(analogue, parms, requireTime = NULL, requireSteps = NULL) {

  lastState <- analogue$state[nrow(analogue), ]
  lastTime <- analogue$time[nrow(analogue)]

  knn <- parms$knnFun(lastState)
  storeIdx <- parms$knnIdxToStoreIdx[knn$idx]

  storeValid <-
    parms$store |>
    dplyr::filter(
      .data$trajId == .data$trajId[storeIdx],
      seq_len(dplyr::n()) >= storeIdx) |>
    dplyr::mutate(time = .data$time - .data$time[1] + lastTime)

  analogue <-
    dplyr::bind_rows(
      analogue,
      storeValid[-1, ]
    )

  timeMax <- max(analogue$time)

  if (!is.null(requireTime) && timeMax >= requireTime) {
    return(analogue)
  }
  if (!is.null(requireSteps) && nrow(analogue) >= requireSteps + 1) {
    return(analogue)
  }

  return(extendAnalogue(analogue, parms, requireTime, requireSteps))
}
