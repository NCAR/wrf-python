; Bug1: This segfaults
l = NewList("fifo")
name = "foo"
ListAppend(l, (/name/))
print(l)
print(l[0])
name = "bar"

; Bug2 Variables disappear
a = addfile("/Users/ladwig/Documents/wrf_files/wrfout_d02_2010-06-13_21:00:00.nc", "r")
b := wrf_user_getvar(a, "slp", -1)
c = NewList("fifo") 
ListAppend(c, (/b/)) 
b := wrf_user_getvar(a, "rh", -1) 
ListAppend(c, (/b/))

print(c[0])
print(c[1]) ; Variables start disappearing
