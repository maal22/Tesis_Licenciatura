import sys
import random
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

##################################################################################################################################
##################################################################################################################################

### CÁLCULO DE DISTANCIAS ENTRE DOS ELEMENTOS DE UNA REJILLA CUADRADA BIDIMENSIONAL DOTADA DE CONDICIONTES DE FRONTERA PERIÓDICAS
def calc_dist(loc1, loc2, longitud):
	distx = loc1[0] - loc2[0]
	disty = loc1[1] - loc2[1]
	dist1 = (distx ** 2 + disty ** 2) ** (1 / 2)
	dist2 = ((distx - longitud) ** 2 + disty ** 2) ** (1 / 2)
	dist3 = ((distx + longitud) ** 2 + disty ** 2) ** (1 / 2)
	dist4 = (distx ** 2 + (disty - longitud) ** 2) ** (1 / 2)
	dist5 = (distx ** 2 + (disty + longitud) ** 2) ** (1 / 2)
	dist6 = ((distx + longitud) ** 2 + (disty + longitud) ** 2) ** (1 / 2)
	dist7 = ((distx - longitud) ** 2 + (disty + longitud) ** 2) ** (1 / 2)
	dist8 = ((distx - longitud) ** 2 + (disty - longitud) ** 2) ** (1 / 2)
	dist9 = ((distx + longitud) ** 2 + (disty - longitud) ** 2) ** (1 / 2)
	lista_distancias = [dist1, dist2, dist3, dist4, dist5, dist6, dist7, dist8, dist9]
	distancia = min(lista_distancias)
	return distancia

### INICIACIÓN DE GRÁFICA PARA SISTEMA CML
def ini_graf(longitud):
	grafica = nx.Graph()
	# Definición de nodos
	num_sitios = longitud ** 2
	for sitio in range(num_sitios):
		grafica.add_node(sitio, u = random.uniform(-1, 1))
	# Definición de bordes
	lista_edges = []
	for sitio1 in range(num_sitios):
		for sitio2 in range(num_sitios):
			if sitio2 == sitio1:
				continue
			coord1 = (sitio1 % longitud, sitio1 // longitud)
			coord2 = (sitio2 % longitud, sitio2 // longitud)
			distancia = calc_dist(coord1, coord2, longitud)
			lista_edges.append((sitio1, sitio2, {'path_distance': distancia}))
	grafica.add_edges_from(lista_edges)
	return grafica

### CÁLCULO DE PRIMEROS VECINOS DE UN NODO EN UNA GRÁFICA CON BORDES DOTADOS DE ATRIBUTOS EQUIVALENTES A LAS DISTANCIAS ENTRE LOS SITIOS QUE LOS DEFINEN
def primeros_vecinos(grafica, nodo):
	# Lista de posibles distancias
	lista_dist = []
	for vec1, datos1 in grafica.adj[nodo].items():
		for keys1, dists1 in datos1.items():
			lista_dist.append(dists1)
	min_dist = min(lista_dist)
	# Lista de sitios a distancia mínima (primeros vecinos)
	lista_1v = []
	for vec2, datos2 in grafica.adj[nodo].items():
		for keys2, dists2 in datos2.items():
			if dists2 == min_dist:
				lista_1v.append(vec2)
	return lista_1v

### MAPEO PHI
def phi(valor_t):
	if -1 <= valor_t < -1 / 3.:
		valor_tt = (-3 * valor_t) - 2
	elif -1 / 3. <= valor_t < 1 / 3.:
		valor_tt = 3 * valor_t
	elif 1 / 3. <= valor_t <= 1:
		valor_tt = (-3 * valor_t) + 2
	return valor_tt

### EVOLUCIÓN TEMPORAL DE LA GRÁFICA
def ev_temp(num_iter, trans, x0, g_acople, grafica, lista_1vecinos, guardar):
	print('INICIO DE EVOLUCION TEMPORAL')
	lista_sitios = list(grafica.nodes())
	tenth_progress = int(num_iter / 10)
	# Guardar evolución de 'u'
	if guardar == True:
		arr_promtemp = np.zeros((num_iter - trans))
		for iteracion1 in range(num_iter):
			if iteracion1 % tenth_progress == 0:
				print('*')
			# Gráfica auxiliar con valores de 'u' de siguiente iteración
			grafica_holder1 = nx.Graph()
			for sitio1a in lista_sitios:
				grafica_holder1.add_node(sitio1a, u = 0)
				sum1vec1 = 0
				for vecino1 in lista_1vecinos[sitio1a][1]:
					dif_u1 = phi(grafica.nodes[vecino1]['u']) - phi(grafica.nodes[sitio1a]['u'])
					sum1vec1 = sum1vec1 + dif_u1
				grafica_holder1.nodes[sitio1a]['u'] = phi(grafica.nodes[sitio1a]['u']) + g_acople * sum1vec1
			# Actualización de gráfica "original"
			for sitio1b in lista_sitios:
				grafica.nodes[sitio1b]['u'] = grafica_holder1.nodes[sitio1b]['u']
			if iteracion1 <= trans - 1:
				continue
			arr_promtemp[iteracion1 - trans] = grafica.nodes[x0]['u']
		print('FIN DE EVOLUCION TEMPORAL')
		return grafica, arr_promtemp
	# Sin guardar evolución de 'u'
	else:
		for iteracion2 in range(num_iter):
			if iteracion2 % tenth_progress == 0:
				print('*')
			# Gráfica auxiliar
			grafica_holder2 = nx.Graph()
			for sitio2a in lista_sitios:
				grafica_holder2.add_node(sitio2a, u = 0)
				sum1vec2 = 0
				for vecino2 in lista_1vecinos[sitio2a][1]:
					dif_u2 = phi(grafica.nodes[vecino2]['u']) - phi(grafica.nodes[sitio2a]['u'])
					sum1vec2 = sum1vec2 + dif_u2
				grafica_holder2.nodes[sitio2a]['u'] = phi(grafica.nodes[sitio2a]['u']) + g_acople * sum1vec2
			# Actualización de gráfica
			for sitio2b in lista_sitios:
				grafica.nodes[sitio2b]['u'] = grafica_holder2.nodes[sitio2b]['u']
		print('FIN DE EVOLUCION TEMPORAL')
		return grafica

##################################################################################################################################
##################################################################################################################################

### DEFINICIÓN DE PARÁMETROS DE SIMULACIÓN
# Selección aleatoria de semilla
s = random.randrange(sys.maxsize)
# Definición de otros parámetros
L = int(input('Ingresa la longitud de la rejilla CML (entero): '))
g = float(input('Ingresa la constante de acoplamiento (real positivo con tres decimales): '))
N_iter = int(input('Ingresa el total de iteraciones (entero): '))
transient = int(input('Ingresa el valor de transient (entero): '))
site_x0 = int(random.randrange(0,int(L**2),1))
N_ensembles = int(input('Ingresa el total de sistemas que conforman el ensamble (entero): '))

### GENERACIÓN DE GRÁFICA Y LISTA CON PRIMEROS VECINOS
# Iniciación de generador de números aleatorios
random.seed(s)
# Definición de gráfica y lista con primeros vecinos
lattice = ini_graf(L)
list_1neighbors = []
for site in range(L**2):
	list1v = primeros_vecinos(lattice, site)
	list_1neighbors.append((site, list1v))
print(lattice.nodes[15]['u'])
print(type(lattice.nodes[15]['u']))

### DISTRIBUCIÓN DE VARIABLE 'U' EN UN SITIO PARTICULAR TRAS MÚLTIPLES ITERACIONES, CONSIDERANDO UN SISTEMA
safe_timeavg = True
lattice, arr_timeavg = ev_temp(N_iter, transient, site_x0, g, lattice, list_1neighbors, safe_timeavg)

fname_timeavg = 'tesis_Egolf_ergodicidad_L%(length)i_g%(coupling).3f_Niter%(iterations).3e_trans%(trans).3e_seed%(seed)i_site%(site)i_timeavg.txt'
dict_fname_timeavg = {'length': L, 'coupling': g, 'iterations': N_iter, 'trans': transient, 'site': site_x0, 'seed': s}
np.savetxt(fname_timeavg % dict_fname_timeavg, arr_timeavg)

print('Evolución de sistema completada')

# Distribución de variable 'u' en un sitio particular tras múltiples iteraciones, considerando un ensamble
safe_ensembleavg = True
arr_ensembleavg = np.zeros((N_ensembles, (N_iter - transient)))
for sys in range(N_ensembles):
	lattice = ini_graf(L)
	lattice, arr_ensemble_holder = ev_temp(N_iter, transient, site_x0, g, lattice, list_1neighbors, safe_ensembleavg)
	arr_ensembleavg[sys] = arr_ensemble_holder
arr_ensembleavg = arr_ensembleavg.flatten()

fname_ensavg = 'tesis_Egolf_ergodicidad_L%(length)i_g%(coupling).3f_Niter%(iterations).3e_trans%(trans).3e_seed%(seed)i_site%(site)i_Nens%(ens)i_ensavg.txt'
dict_fname_ensavg = {'length': L, 'coupling': g, 'iterations': N_iter, 'trans': transient, 'site': site_x0, 'ens': N_ensembles, 'seed': s}
np.savetxt(fname_ensavg % dict_fname_ensavg, arr_ensembleavg)

print('Evolución de ensamble completada')

# Iniciación de rejilla CML para prueba de self-averaging
L2 = int(2*L)
lattice2 = ini_graf(L2)
list_1neighbors2 = []
for site2 in range(L2**2):
	list1v2 = primeros_vecinos(lattice2, site2)
	list_1neighbors2.append((site2, list1v2))

# Distribución de variable 'u' en un sistema a un tiempo fijo
safe_selfavg = False
lattice2 = ev_temp(N_iter, transient, site_x0, g, lattice2, list_1neighbors2, safe_selfavg)
arr_selfavg = np.zeros(L2**2)
for site in range(L2**2):
	arr_selfavg[site] = lattice2.nodes[site]['u']

fname_selfavg = 'tesis_Egolf_ergodicidad_L%(length)i_g%(coupling).3f_Niter%(iterations).3e_trans%(trans).3e_seed%(seed)i_selfavg.txt'
dict_fname_selfavg = {'length': L, 'coupling': g, 'iterations': N_iter, 'trans': transient, 'seed': s}
np.savetxt(fname_selfavg % dict_fname_selfavg, arr_selfavg)

print('Prueba de self-averaging completada')

# Gráfica con resultados de ergodicidad y self-averaging
plt.figure(1)
fig1, (ax1A, ax1B, ax1C) = plt.subplots(nrows = 1, ncols = 3, figsize = (30,10))
plt.tight_layout(pad=4, h_pad=4, w_pad=6)

hist_time, bins_time = np.histogram(arr_timeavg, range = (-1, 1))
hist_ens, bins_ens = np.histogram(arr_ensembleavg, range = (-1, 1))
hist_self, bins_self = np.histogram(arr_selfavg, range = (-1,1))

ax1A.hist(bins_time[:-1], bins_time, weights = hist_time, density = True)
ax1A.set_title('Distribución de ' + r'$u_{\vec{x}_{0}}^{t}$' + ' con un sistema de %(lo)i x %(lo)i sitios ' % {'lo': L} + r'$(\vec{x}_{0} = %(siteref)i) $' % {'siteref': site_x0}, size=16)
ax1A.set_ylabel('dP(u)', size=15)
ax1A.set_xlabel('u', size=15)
for tick in ax1A.xaxis.get_major_ticks():
	tick.label.set_fontsize(14)
for tick in ax1A.yaxis.get_major_ticks():
	tick.label.set_fontsize(14)

ax1B.hist(bins_ens[:-1], bins_ens, weights = hist_ens, density = True)
ax1B.set_title('Distribución de ' + r'$u_{\vec{x}_{0}}^{t}$' + ' con un ensamble de %(Nens)i sistemas de %(lo)i x %(lo)i sitios' % {'Nens': N_ensembles, 'lo': L}, size=16)
ax1B.set_ylabel('dP(u)', size=15)
ax1B.set_xlabel('u', size=15)
for tick in ax1B.xaxis.get_major_ticks():
	tick.label.set_fontsize(14)
for tick in ax1B.yaxis.get_major_ticks():
	tick.label.set_fontsize(14)

ax1C.hist(bins_self[:-1], bins_self, weights = hist_self, density = True)
ax1C.set_title('Distribución de ' + r'$u_{\vec{x}}^{t_{0}}$' + ' sobre un sistema de %(lo)i x %(lo)i sitios ' % {'lo': L2} + r'$(t_{0} = %(tiempo).2e) $' % {'tiempo': N_iter}, size=16)
ax1C.set_ylabel('dP(u)', size=15)
ax1C.set_xlabel('u', size=15)
for tick in ax1C.xaxis.get_major_ticks():
	tick.label.set_fontsize(14)
for tick in ax1C.yaxis.get_major_ticks():
	tick.label.set_fontsize(14)

imgname = 'tesis_Egolf_ergodicidad_L%(length)i_g%(coupling).3f_Niter%(iterations).3e_trans%(trans).3e_seed%(seed)i_site%(site)i_Nens%(ens)i_Graph.png'
dict_imgname = {'length': L, 'coupling': g, 'iterations': N_iter, 'trans': transient, 'site': site_x0, 'ens': N_ensembles, 'seed': s}
plt.savefig(imgname % dict_imgname)

print('Programa concluido')

