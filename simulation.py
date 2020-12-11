# LEPL1501 annee 2020 
# Simualation: oscillation d'une maquette de grue flottante en fonction du temps
# Groupe 11.11
# Alexis Airson, Ilias Asqiriba, Sacha Aoust, Emilie Art, Lou-Anne Bertrand, Cedric Amerijckx
# Code inspire par Charles Pecheur 2018

#Par respect des conventions, les unites sont les metres, les kilogrammes, les secondes, les radians et les joules

#L'origine du repere se trouve au milieu de la longueur de la plateforme et au niveau de la ligne de flottaison 

from math import *
import matplotlib.pyplot as plot
import numpy as np

#Constantes & parametres:

g = 9.81     # constante gravitationnelle
rho_eau = 1000   # masse volumique de l'eau

masse_plateforme = 7.5
masse_grue = 1.5
masse_charge = 0.2 # 200g (voir instructions)
somme_masses  = masse_plateforme + masse_grue + masse_charge

h1 = 0.1 # hauteur de la plateforme 
l = 0.6  # longueur et largeur de la plateforme carree
distance = 0.7 # 40cm comme demande entre la plateforme et la structure + la moitie de la longueur soit 30 cm
distance_max = 0.8 # on prevoit d'aller un poil plus loin 

I = 1.18 # inertie de l'engin (obtenue avec fusion)
D = 1 # coefficient d'amortissement de l'eau

#Parametres de la simulation :

dt = 0.001        # pas (dt) 
end = 10          # duree

t = np.arange(0, end, dt) # Tableau principal, stocke le temps en secondes
theta = np.empty_like(t) # Tableau de l'angle, en radians
omega = np.empty_like(t) # Tableau de la vitesse angulaire, en rad/s
a = np.empty_like(t) # Tableau de l'acceleration angulaire, en rad/s**2
d = np.empty_like(t) # Tableau de la distance, en metres
c = np.empty_like(t) # Tableau des couples, en Newton-metre
Ek = np.empty_like(t) # Tableau de l'energie cinetique, en Joules
Eg = np.empty_like(t) # Tableau de l'energie potentielle, en Joules
Em = np.empty_like(t) # Tableau de l'energie mecanique, en Joules
theta[0] = 0 
omega[0] = 0 #conditions initiales 

# calcul de hc (sachant que la poussee d'archimede compense toutes les forces poids):

def compute_hc():
    
    ''' pre : None
        post : la fonction renvoie la hauteur hc qui represente la hauteur immergee'''
    
    hc = (g*somme_masses)/(g*l**2*rho_eau)
    return(hc)

# calcul du centre de gravite initial (cfr livret physique S1) :

def cg_initial():
    
    ''' pre : None
        post : la fonction renvoie les coordonees du centre gravite initial avec les variables xg0 et zg0 '''
    
    xg0 = -4.866/100
    zg0 = 5.751/100  # obtenus via Fusion 360
    return(xg0,zg0)

# calcul du nouveau centre de poussee :

def centre_poussee(theta,hc) :
    
    ''' pre : un angle theta (theta[i] dans la simulation) et hc (obtenu via la fonction compute)
        post : la fonction renvoie les coordonnees du centre de poussee pour un certain angle theta'''
    
    xa = l/2
    za = -hc     # pour comprendre ceci, lire la partie physique 
    
    h1= hc + tan(theta)*l/2
    h2= hc - tan(theta)*l/2
    
    l_C=(l*(h1+2*h2))/(3*(h1+h2))
    H_C=(h1**2+h1*h2+h2**2)/(3*(h1+h2))
    
    xc = xa - l_C
    zc = za + H_C
    
    XC = xc*cos(theta) + zc*sin(theta)
    ZC = -xc*sin(theta) + zc*cos(theta)
    
    return(XC,ZC)

# calcul de l'angle theta max :

def theta_max(hc):
    
    ''' pre : hc (obtenu via la fonction compute)
        post : la fonction renvoie les angles de soulevement et de submersion'''
    
    return (atan((h1-hc)/(l/2)) , atan(hc/(l/2))) #cfr exercices de physique 

# distance qui evolue au fil du temps :

def distance_temps():
    
    ''' pre : None
        post : la fonction renvoie la distance, correspondant a un certain pourcentage de la distance maximale (cas d'une distance variable)'''
    
    for x in range(0,len(t)-1):
        d[x]=x/(len(t))*dmax  # La distance correspond a un certain pourcentage de la distance maximale
        
# Expression de l'angle predit :

def angle():
    
    ''' pre : None
        post : la fonction affiche l'angle predit en degres et en radians sur la console'''
    
    last = end/dt - 1 # ceci represente le dernier angle pris en compte dans la simulation
    
    print(str(abs(theta[int(last)])) + " radians") # affiche l'angle predit sur la console
    
    last = abs(theta[int(last)]) * 180 / pi # convertit l'angle en degres
    
    print(str(last) + ' degres') # affiche l'angle predit sur la console 

# recuperation variables :

def parametres():
    
    ''' pre : None
        post : la fonction fait appel aux fonctions precedentes et renvoie tous les parametres calcules dans ces fonctions'''
    
    hc = compute_hc()
    cg_x = cg_initial()[0]
    cg_z = cg_initial()[1]
    cp_x = centre_poussee(0,hc)[0]
    cp_z = centre_poussee(0,hc)[1]
    theta_submersion,theta_soulevement = theta_max(hc)
    
    return (hc,cg_x,cg_z,cp_x,cp_z,theta_submersion,theta_soulevement)

# definition de la simulation avec une distance fixe (distance_max) :

def simulation(somme_masses):
    
    ''' pre : la valeur correspondant a la masse totale du systeme (evite la creation d'une variable locale)
        post : la fonction effectue la simulation et remplit les tableaux definis plus haut'''
    
    #conditions initiales
    
    xg0 = cg_initial()[0] 
    zg0 = cg_initial()[1]
    
    for i in range(len(t)-1):
        
        d[i] = distance_max # distance fixe 
        
        cp_x = centre_poussee(theta[i],hc)[0] # calcul du centre de poussee 
        
        xg = xg0*cos(theta[i]) + sin(theta[i])*zg0 #on fait simplement subir une rotation au cg
        zg = -xg0*sin(theta[i]) + cos(theta[i])*zg0
        
        Eg[i] = somme_masses*g*zg # complete le tableau de l'energie potentielle 
        
        Ca = masse_charge*g*d[i] # calcul du couple lie a la charge 
        Cr = g*(masse_plateforme + masse_grue)*(cp_x - xg) # calcul du couple de redresssement 
        C = Ca - Cr # calcul du couple global (-Cr pcq les couples ne sont pas dans le meme sens)
        
        Ek[i] = I*(omega[i])**2*(1/2) # complete le tableau de l'energie cinetique 
        
        a[i]= (-D*omega[i]+C)/I # complete le tableau de l'acceleration angulaire 
        omega[i+1] = omega[i] + a[i]*dt # complete le tableau de la vitesse angulaire 
        theta[i+1] = theta[i] + omega[i]*dt # complete le tableau de l'angle d'inclinaison 
        #cfr l'equation differentielle
        
        Em[i] = Ek[i] + Eg[i] # complete le tableau de l'energie mecanique 
        
        last = end/dt - 1 # ceci represente le dernier angle pris en compte dans la simulation
        
        Eg[int(last)] = zg*g*somme_masses # complete le tableau de l'energie potentielle
    
        Em[int(last)] = Eg[int(last)] + Ek[int(last)] # complete le tableau de l'energie mecanique 
        
# definition de la simulation avec une distance evoluant au fil du temps :

def simulation_distance_variable(somme_masses):
    
    #conditions initiales
    
    xg0 = cg_initial()[0] 
    zg0 = cg_initial()[1]
    
    for i in range(len(t)-1):
        
        d[i]=i/(len(t)-1)*distance_max
        
        cp_x = centre_poussee(theta[i],hc)[0] # calcul du centre de poussee 
        
        xg = xg0*cos(theta[i]) + sin(theta[i])*zg0 #on fait simplement subir une rotation au cg
        zg = -xg0*sin(theta[i]) + cos(theta[i])*zg0
        
        Eg[i] = somme_masses*g*zg # complete le tableau de l'energie potentielle 
        
        Ca = masse_charge*g*d[i] # calcul du couple lie a la charge 
        Cr = g*(masse_plateforme + masse_grue)*(cp_x - xg) # calcul du couple de redresssement 
        C = Ca - Cr # calcul du couple global (-Cr pcq les couples ne sont pas dans le meme sens)
        
        Ek[i] = I*(omega[i])**2*(1/2) # complete le tableau de l'energie cinetique 
        
        a[i]= (-D*omega[i]+C)/I # complete le tableau de l'acceleration angulaire 
        omega[i+1] = omega[i] + a[i]*dt # complete le tableau de la vitesse angulaire 
        theta[i+1] = theta[i] + omega[i]*dt # complete le tableau de l'angle d'inclinaison 
        #cfr l'equation differentielle
        
        Em[i] = Ek[i] + Eg[i] # complete le tableau de l'energie mecanique 
        
        last = end/dt - 1 # ceci represente le dernier angle pris en compte dans la simulation
        
        Eg[int(last)] = zg*g*somme_masses # complete le tableau de l'energie potentielle
    
        Em[int(last)] = Eg[int(last)] + Ek[int(last)] # complete le tableau de l'energie mecanique
        
# definition des graphiques

# angle, vitesse et acceleration par rapport au temps (distance fixe)

def graphiques():
    
    ''' pre : None
        post : Renvoie les graphiques de l'acceleration angulaire, la vitesse angulaire et l'angle d'inclinaison pour une distance fixe'''
    
    plot.suptitle('angle, vitesse et acceleration par rapport au temps (avec distance fixe)')
    plot.subplot(3,1,1)
    plot.ylabel(r'angle [rad]')
    plot.plot(t,theta, label="angle theta",color='blue')
    plot.hlines([theta_submersion,-theta_submersion],0,end,color='red',label='submersion')
    plot.hlines([theta_soulevement,-theta_soulevement],0,end,color='green',label='soulevement')
    plot.legend(loc = "upper right", prop={'size': 8})
    
    plot.subplot(3,1,2)
    plot.ylabel(r'ω [rad/s]')
    plot.plot(t,omega, label="v angulaire omega",color='black')
    plot.legend()
    
    plot.subplot(3,1,3)
    plot.ylabel(r"a [rad/s^2]")
    plot.xlabel(r'time [s]')
    plot.plot(t,a, label="a angulaire",color='purple')
    plot.legend()
    plot.show()    
    
    plot.title('theta en fonction du temps')
    plot.ylabel(r'angle [rad]')
    plot.xlabel(r'time [s]')
    plot.plot(t,theta, label="angle theta",color='blue')   
    plot.hlines([theta_submersion,-theta_submersion],0,end,color='red',label='submersion')
    plot.hlines([theta_soulevement,-theta_soulevement],0,end,color='green',label='soulevement')
    plot.legend()
    plot.autoscale()
    plot.show()
   
def graphiques_distance_variable():
    
    ''' pre : None
        post : Renvoie les graphiques de l'acceleration angulaire, la vitesse angulaire et l'angle d'inclinaison pour une distance variable'''
    
    plot.suptitle('angle, vitesse et acceleration par rapport au temps (avec distance variable)')
    plot.subplot(3,1,1)
    plot.ylabel(r'angle [rad]')
    plot.plot(t,theta, label="angle theta",color='blue')
    plot.hlines([theta_submersion,-theta_submersion],0,end,color='red',label='submersion')
    plot.hlines([theta_soulevement,-theta_soulevement],0,end,color='green',label='soulevement')
    plot.legend(loc = "upper right", prop={'size': 8})
    
    plot.subplot(3,1,2)
    plot.ylabel(r'ω [rad/s]')
    plot.plot(t,omega, label="v angulaire omega",color='black')
    plot.legend()
    plot.subplot(3,1,3)
    plot.ylabel(r"a [rad/s^2]")
    plot.xlabel(r'time [s]')
    plot.plot(t,a, label="a angulaire",color='purple')
    plot.legend()
    plot.show()
      
    plot.title('theta en fonction du temps')
    plot.ylabel(r'angle [rad]')
    plot.xlabel(r'time [s]')
    plot.plot(t,theta, label="angle theta",color='blue') 
    plot.hlines([theta_submersion,-theta_submersion],0,end,color='red',label='submersion')
    plot.hlines([theta_soulevement,-theta_soulevement],0,end,color='green',label='soulevement')
    plot.legend()
    plot.autoscale()
    plot.show()
    
def graphique_Ek():
    
    ''' pre : None
        post : renvoie le graphique de l'energie cinetique en fonction du temps '''
    
    plot.title('energie cinetique en fonction du temps')
    plot.ylabel(r'energy [J]')
    plot.xlabel(r'time [s]')
    plot.plot(t,Ek,color = 'orange')
    plot.show()

def graphique_Eg():
    
    ''' pre : None
        post : renvoie le graphique de l'energie potentielle en fonction du temps '''
    
    plot.title('energie potentielle en fonction du temps')
    plot.ylabel(r'energy [J]')
    plot.xlabel(r'time [s]')
    plot.plot(t,Eg,color = 'brown')
    plot.show()
    
def graphique_Em():
    
    ''' pre : None
        post : renvoie le graphique de l'energie mecanique en fonction du temps '''
    
    plot.title('energie mecanique en fonction du temps')
    plot.ylabel(r'energy [J]')
    plot.xlabel(r'time [s]')
    plot.plot(t,Em,color = 'cyan')
    plot.show()

def diagramme():
    
    ''' pre : None
        post : renvoie le graphique du diagramme de phase (vitesse angulaire en fonction de l'angle d'inclinaison) '''
    
    plot.title('diagramme de phase')
    plot.ylabel(r'omega[rad/s]')
    plot.xlabel(r'theta[rad]')
    plot.plot(theta,omega,color='magenta')
    plot.show()
    
#Programme principal

if __name__=='__main__': # verifie que le code est bien effectue depuis ce fichier 
    n=0
    while n < 1 : # Possibilite d'effectuer plusieurs fois la simulation
        theta[0] = 0
        omega[0] = 0
        hc,cg_x,cg_z,cp_x,cp_z,theta_submersion,theta_soulevement = parametres() # initialisation des parametres 
        simulation(somme_masses)
        graphiques()
        angle()
        graphique_Ek()
        graphique_Eg()
        graphique_Em()
        diagramme()
        hc,cg_x,cg_z,cp_x,cp_z,theta_submersion,theta_soulevement = parametres() #reinitialisation des parametres 
        simulation_distance_variable(somme_masses)
        graphiques_distance_variable()
        angle()
        graphique_Ek()
        graphique_Eg()
        graphique_Em()
        diagramme()
        n += 1
