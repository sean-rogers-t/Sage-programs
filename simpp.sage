def mincp(n,k,m):
	W=WeylGroup('C'+str(n), prefix="s")
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
	Wpm=[]
	for i in range(len(Wp)):
		if Wp[i].length() == m:
			Wpm.append(Wp[i])
	return Wpm
	
	
def mincp_total(n,k):
	W=WeylGroup('C'+str(n), prefix="s")
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

def cp_long(n,k):
	cp_length=[]
	for i in range(len(mincp_total(n,k))):
		cp_length.append(mincp_total(n,k)[i].length())
	ml=max(cp_length)
	ml_index=cp_length.index(ml)
	cpl=mincp_total(n,k)[ml_index]
	return cpl
	
def long(n):
	W=WeylGroup('C'+str(n), prefix="s")
	simp=W.simple_reflections()
	simp1=list(simp)
	x=W.long_element()
	return x
def PD(n,k,e):
	w_nk=cp_long(n,k)
	w0=long(n)
	return w0*e*w_nk*w0
def sr(n,elem):
	W=WeylGroup('C'+str(n), prefix="s")
	WA=WeylGroup('A'+str((2*n)-1),prefix='r')
	sc=W.simple_reflections()
	sa=WA.simple_reflections()
	str_elem=str(elem)
	old_elem=[]
	for i in range(len(str_elem)):
		if i%3==1:
			old_elem.append(int(str_elem[i]))
	new_elem=[]
	for i in range(len(old_elem)):
		if old_elem[i]!=n:
			new_elem.append(old_elem[i])
			new_elem.append((2*n)-old_elem[i])
		if old_elem[i]==n:
			new_elem.append(n)
	ref=WA[0]
	for i in range(len(new_elem)):
		ref=ref*sa[new_elem[i]]
	return ref