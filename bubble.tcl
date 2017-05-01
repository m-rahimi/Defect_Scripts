wrapmode cell

foreach { x0 y0 z0 } $bubbleCenter { break }

proc calcforces {step unique K R1 R2 feq} {

  global x0 y0 z0 
  set Pi 3.14159265359

  if { $step % $feq == 0 } { cleardrops } 

  while {[nextatom]} { 
  
    set rvec [getcoord]
    foreach { x y z } $rvec { break }
    set RP2 [expr {(($x-$x0)*($x-$x0) + ($y-$y0)*($y-$y0))}]
    set R  [expr {sqrt($RP2)}]
         
    if { $R < $R1 } {
        
      set ff     [expr {(($K * $Pi / $R1) * sin($Pi * $R / $R1))}]
      set forceX [expr {($ff * (($x-$x0) / $R))}]
      set forceY [expr {($ff * (($y-$y0) / $R))}]
      set forceZ 0.000
      addforce "$forceX $forceY $forceZ"

    } elseif { $R > $R2 } {
        
      dropatom

    }
  }
}
