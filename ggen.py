#!/usr/bin/python
from sys import argv
from random import randint

n = int(argv[1])
k = int(argv[2])

print n+2, n
for i in xrange(n):
	print "i%d" % (i+1)

l = [[] for i in xrange(2+n)]
l[0].append((1,100))
l[1].append((0,100))
conn=set()
conn.add((0,1))
conn.add((1,0))
for i in xrange(k):
	a,b = randint(0,n+1),randint(0,n+1)
	while a==b or (a,b) in conn:
		b = randint(0,n+1)
	d = randint(1, 100)
	l[a].append((b,d))
	l[b].append((a,d))
	conn.add((a,b))
	conn.add((b,a))

for i in l:
	print len(i), 0, 0
	for j,k in sorted(i):
		print j, k,
	print
