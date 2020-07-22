def elem(k,h):
	s=range(1,k+1)
	power_set=[[]]
	power_set_h=[]
	for elem in s:
		#iterate overs ubsets so far
		for sub_set in power_set:
			#add a new subset consiting of the subset at hand added elem
			power_set=power_set+[list(sub_set)+[elem]]
	for i in range(len(power_set)):
		if len(power_set[i])==h:
			power_set_h.append(power_set[i])
	return power_set_h

def chev_tot_d(n,k,h):
	R=PolynomialRing(ZZ,['x%s'%p for p in range(1,n+1)]) #gives polynomial ring in n variables
	R.inject_variables() #makes the variables useable
	Gen=R.gens() #gives list of variables
	lst=[(Gen[-i-1]-Gen[-i-2])^2 for i in range(1,k-1)]
	lst.append((Gen[-2]+Gen[-1]-Gen[-3])^2) #gives image under associated Springer map
	lst.append((Gen[-1]-Gen[-2])^2)

	full_elem=[] #expands h elemenatary polynomial in Springer variables from lst
	for i in range(len(elem(k,h))):
		f=1
		for j in range(len(elem(k,h)[i])):
			f=f*lst[elem(k,h)[i][j]-1]
		full_elem.append(f)
	g=0
	for i in range(len(full_elem)):
		#adds all elements in full_elem
		g=g+full_elem[i]
	
	co=g.coefficients() #gives list of coefficients from g
	exp=g.exponents()   #gives list of exponents as tuples from g

	chev_exp=[] #puts exponents in proper form for finchev function
	for i in range(len(exp)):
		new_exp=[]
		for j in range(len(exp[0])):
			for k in range(exp[i][j]):
				new_exp.append(j+1)
		chev_exp.append(new_exp)

	chev_total=[] #gives list of coefficients with simple reflection numbers
	for i in range(len(co)):
		chev_total.append([co[i],chev_exp[i]])
	
	return chev_total

def special_element(n,k):
	#gives image of "x_(n-k+1)*...*x_n"
	R=PolynomialRing(ZZ,['x%s'%p for p in range(1,n+1)]) #gives polynomial ring in n variables
	R.inject_variables() #makes the variables useable
	Gen=R.gens() #gives list of variables
	lst=[(Gen[-i-1]-Gen[-i-2]) for i in range(2,k)]
	lst.append((Gen[-2]+Gen[-1]-Gen[-3])) #gives image under associated Springer map
	lst.append((Gen[-1]-Gen[-2]))

	
	g=1
	for i in range(len(lst)):
		#product of all elements in lst
		g=g*lst[i]
	
	co=g.coefficients() #gives list of coefficients from g
	exp=g.exponents()   #gives list of exponents as tuples from g

	chev_exp=[] #puts exponents in proper form for finchev function
	for i in range(len(exp)):
		new_exp=[]
		for j in range(len(exp[0])):
			for k in range(exp[i][j]):
				new_exp.append(j+1)
		chev_exp.append(new_exp)

	chev_total=[] #gives list of coefficients with simple reflection numbers
	for i in range(len(co)):
		chev_total.append([co[i],chev_exp[i]])
	
	return chev_total
def boom_d(n,k,h):
	chevy=chev_tot_d(n,k,h)
	deco=[chevy[i][1] for i in range(len(chevy))]
	cof=[chevy[i][0] for i in range(len(chevy))]

	multy_chev=[multlist(cof[i],finchev_d(n,deco[i])) for i in range(len(chevy))]
	one_list=[]
	for i in range(len(multy_chev)):
		for j in range(len(multy_chev[i])):
			one_list.append(multy_chev[i][j])
	return one_list
def boom_spec(n,k):
	chevy=special_element(n,k)
	deco=[chevy[i][1] for i in range(len(chevy))]
	cof=[chevy[i][0] for i in range(len(chevy))]

	multy_chev=[multlist(cof[i],finchev_d(n,deco[i])) for i in range(len(chevy))]
	one_list=[]
	for i in range(len(multy_chev)):
		for j in range(len(multy_chev[i])):
			one_list.append(multy_chev[i][j])
	return one_list
def cohom_spec(n,k):
	pre_add=boom_spec(n,k)
	pre_set=[pre_add[i][1] for i in range(len(pre_add))]
	cox=list(set(pre_set))
	tot_list=[]
	for i in range(len(cox)):
		cof_i=0
		for j in range(len(pre_set)):
			if cox[i]==pre_set[j]:
				cof_i=cof_i+pre_add[j][0]
		tot_list.append([cof_i,cox[i]])
	return tot_list

def cohom_0_spec(n,k):
	non_0=cohom_spec(n,k)
	coh_0=[]
	for i in range(len(non_0)):
		if non_0[i][0]!=0:
			coh_0.append(non_0[i])
	return coh_0
def cohom_d(n,k,h):
	pre_add=boom_d(n,k,h)
	pre_set=[pre_add[i][1] for i in range(len(pre_add))]
	cox=list(set(pre_set))
	tot_list=[]
	for i in range(len(cox)):
		cof_i=0
		for j in range(len(pre_set)):
			if cox[i]==pre_set[j]:
				cof_i=cof_i+pre_add[j][0]
		tot_list.append([cof_i,cox[i]])
	return tot_list
	

def cohom_0_d(n,k,h):
	non_0=cohom_d(n,k,h)
	coh_0=[]
	for i in range(len(non_0)):
		if non_0[i][0]!=0:
			coh_0.append(non_0[i])
	return coh_0
		
	