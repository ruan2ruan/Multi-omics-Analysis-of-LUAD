display.progress = function ( index, totalN, breakN=20) {

if ( index %% ceiling(totalN/breakN)  ==0  ) {
      cat(paste(round(index*100/totalN), "% ", sep=""))
}

}    