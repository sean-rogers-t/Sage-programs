def mincp_a(n,k,m):
	W=WeylGroup('A'+str(n-1), prefix='r')
	simp=W.simple_reflections()
	simp1=list(simp)
	#takes simple reflections out of dictionary
	simpp=[]
	for i in range(len(simp1)):
		if i!=(n-(k+1)):
			simpp.append(simp1[i])
	#gets rid of important simple reflection
	Wp=[]
	#will be minimal length coset representatives
	for i in range(len(W)):
		if all((W[i]*simpp[j]).length()>W[i].length() for j in range(len(simpp))) == True:
			Wp.append(W[i])
	Wpm=[]
	#will give minimal length coset representatives of the length m
	for i in range(len(Wp)):
		if Wp[i].length()==m:
			Wpm.append(Wp[i])
	return Wpm

def mincp_total_a(n,k):
	W=WeylGroup('A'+str(n-1), prefix="r")
	simp=W.simple_reflections()
	simp1=list(simp)
	simpp=[]
	for i in range(len(simp1)):
		if i !=(n-(k+1)):
			simpp.append(simp1[i])
	Wp=[]
	for i in range(len(W)):
		if all((W[i]*simpp[j]).length()>W[i].length() for j in range(len(simpp))) == True:
			Wp.append(W[i])
	return Wp

def cp_long_a(n,k):
	cp_length=[]
	for i in range(len(mincp_total_a(n,k))):
		cp_length.append(mincp_total_a(n,k)[i].length())
	ml=max(cp_length)
	ml_index=cp_length.index(ml)
	cpl=mincp_total(n,k)[ml_index]
	return cpl
	
def long_a(n):
	W=WeylGroup('A'+str(n-1), prefix="r")
	simp=W.simple_reflections()
	simp1=list(simp)
	x=W.long_element()
	return x
def PD_a(n,k,e):
	w_nk=cp_long_a(n,k)
	w0=long_a(n)
	return w0*e*w_nk*w0
	
	