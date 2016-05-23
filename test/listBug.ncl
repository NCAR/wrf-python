l = NewList("fifo")
name = "foo"
ListAppend(l, (/name/))
print(l)
print(l[0])
name = "bar"
