import vacuumms

c = vacuumms.Configuration('')

o = vacuumms.DDX(c)
o.execute() 

cavs = o.getResult()

dist = vacuumms.csd(cavs)





