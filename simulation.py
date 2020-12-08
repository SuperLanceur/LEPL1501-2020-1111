# LEPL1501 annee 2020 
# Simualation: oscillation d'une maquette de grue flottante en fonction du temps
# Groupe 11.11
# Alexis Airson, Ilias Asqiriba, Sacha Aoust, Emilie Art, Lou-Anne Bertrand, Cedric Amerijckx
# Code inspire par Charles Pecheur 2018

#Par respect des conventions, les unites sont les metres, les kilogrammes, les secondes et les radians 

#L'origine du repere se trouve au milieu de la longueur de la plateforme et au niveau de la ligne de flottaison 

from math import *
import matplotlib.pyplot as plot
import numpy as np

#Constantes & parametres:

g = 9.81     # constante gravitationnelle
rho_eau = 1000   # masse volumique de l'eau

masse_plateforme = 6
masse_grue = 1
masse_charge = 0.2 # 200g (voir instructions)
somme_masses  = masse_plateforme + masse_grue + masse_charge

h1 = 0.1 # hauteur de la plateforme 
l = 0.6  # longueur et largeur de la plateforme carree
distance = 0.7 # 40cm comme demandé entre la plateforme et la structure + la moitié de la longueur soit 30 cm
distance_max = 0.8 # on prevoit d'aller un poil plus loin 

I = 0.74 # inertie de l'engin (obtenue avec fusion)
D = 1 # coefficient d'amortissement de l'eau

#Paramètres de la simulation :

dt = 0.001        # pas (dt) 
end = 10          # duree

t = np.arange(0, end, dt) # Tableau principal, stocke le temps
theta = np.empty_like(t) # Tableau de l'angle, en rad
omega = np.empty_like(t) # Tableau de la vitesse angulaire, en rad/s
a = np.empty_like(t) # Tableau de l'accélération angulaire, en rad/s**2
d = np.empty_like(t) # Tableau de la distance, en mètres
c = np.empty_like(t) # Tableau des couples, en Newton-mètre
theta[0] = 0
omega[0] = 0

# calcul de hc (sachant que la poussee d'archimede compense toutes les forces poids):

def compute_hc():
    hc = (g*somme_masses)/(g*l**2*rho_eau)
    return(hc)

# calcul du centre de gravite initial (cfr livret physique S1) :

def cg_initial():
    xg0 = -4.866/100
    zg0 = 5.751/100
    return(xg0,zg0)

# calcul du nouveau centre de poussée :

def centre_poussee(theta,hc) :
    
    xa = l/2
    za = -hc     # pour comprendre ceci, lire README 1
    
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
    return (atan((h1-hc)/(l/2)) , atan(hc/(l/2))) #cfr exercices de physique 

# distance qui evolue au fil du temps :

def distance_temps():
    for x in range(0,len(t)-1):
        d[x]=x/(len(t))*dmax

# recuperation variables :

def parametres():
    hc = compute_hc()
    cg_x = cg_initial()[0]
    cg_z = cg_initial()[1]
    cp_x = centre_poussee(0,hc)[0]
    cp_z = centre_poussee(0,hc)[1]
    theta_submersion,theta_soulevement = theta_max(hc)
    
    return (hc,cg_x,cg_z,cp_x,cp_z,theta_submersion,theta_soulevement)

# definition de la simulation avec une distance fixe (distance_max) :

def simulation():
    
    #conditions initiales
    
    xg0 = cg_initial()[0] 
    zg0 = cg_initial()[1]
    
    for i in range(len(t)-1):
        
        d[i] = distance_max
        
        cp_x = centre_poussee(theta[i],hc)[0]
        xg = xg0*cos(theta[i]) + sin(theta[i])*zg0 #on fait simplement subir une rotation au cg 
        
        Ca = masse_charge*g*d[i]
        Cr = g*(masse_plateforme + masse_grue)*(cp_x - xg)
        C = Ca - Cr
        
        a[i]= (-D*omega[i]+C)/I #cfr l'equation differentielle
        omega[i+1] = omega[i] + a[i]*dt
        theta[i+1] = theta[i] + omega[i]*dt
        

# definition de la simulation avec une distance évoluant au fil du temps :

def simulation_distance_variable():
    
    #conditions initiales
    
    xg0 = cg_initial()[0] 
    zg0 = cg_initial()[1]
    
    for i in range(len(t)-1):
        
        d[i]=i/(len(t)-1)*distance_max
        
        cp_x = centre_poussee(theta[i],hc)[0]
        xg = xg0*cos(theta[i]) + sin(theta[i])*zg0 #on fait simplement subir une rotation au cg 
        
        Ca = masse_charge*g*d[i]
        Cr = g*(masse_plateforme + masse_grue)*(cp_x - xg)
        C = Ca - Cr
        
        a[i]= (-D*omega[i]+C)/I #cfr l'equation differentielle 
        omega[i+1] = omega[i] + a[i]*dt
        theta[i+1] = theta[i] + omega[i]*dt
# definition des graphiques

# angle, vitesse et accélération par rapport au temps (distance fixe)

def graphiques():
    
    plot.suptitle('angle, vitesse et accélération par rapport au temps (avec distance fixe)')
    plot.subplot(3,1,1)
    plot.ylabel(r'angle [rad]')
    plot.plot(t,theta, label="angle theta",color='blue')
    plot.hlines([theta_submersion,-theta_submersion],0,end,color='red',label='submersion')
    plot.hlines([theta_soulevement,-theta_soulevement],0,end,color='green',label='soulèvement')
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
    plot.hlines([theta_soulevement,-theta_soulevement],0,end,color='green',label='soulèvement')
    plot.legend()
    plot.autoscale()
    theta_end=theta[len(t)-1]
    plot.show()
   
def graphiques_distance_variable():
    
    plot.suptitle('angle, vitesse et accélération par rapport au temps (avec distance variable)')
    plot.subplot(3,1,1)
    plot.ylabel(r'angle [rad]')
    plot.plot(t,theta, label="angle theta",color='blue')
    plot.hlines([theta_submersion,-theta_submersion],0,end,color='red',label='submersion')
    plot.hlines([theta_soulevement,-theta_soulevement],0,end,color='green',label='soulèvement')
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
    plot.hlines([theta_soulevement,-theta_soulevement],0,end,color='green',label='soulèvement')
    plot.legend()
    plot.autoscale()
    theta_end=theta[len(t)-1]
    plot.show()
    
#Programme principal
    
n=0
while n < 1 : #Possibilité d'effectuer plusieurs fois la simulation
    theta[0] = 0
    omega[0] = 0
    hc,cg_x,cg_z,cp_x,cp_z,theta_submersion,theta_soulevement = parametres()
    simulation()
    graphiques()
    hc,cg_x,cg_z,cp_x,cp_z,theta_submersion,theta_soulevement = parametres()
    simulation_distance_variable()
    graphiques_distance_variable()
    n += 1
