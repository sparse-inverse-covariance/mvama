netplot <-
function(g,layout=layout.kamada.kawai,seed=123,...)
  {
  #... could be vertex.label.color="black" or edge.color="black" ....
  set.seed(seed)
  plot(g,layout=layout.kamada.kawai,...)
  }

