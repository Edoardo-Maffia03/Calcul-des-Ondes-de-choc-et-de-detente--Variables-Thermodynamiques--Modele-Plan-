# -*- coding: utf-8 -*-
"""
===============================================================================
   ____            __                 __              
  / __ \____  ____/ /__  _____   ____/ /__      _____/ /_  ____  _____
 / / / / __ \/ __  / _ \/ ___/  / __  / _ \    / ___/ __ \/ __ \/ ___/
/ /_/ / / / / /_/ /  __(__  )  / /_/ /  __/   / /__/ / / / /_/ / /__  
\____/_/_/_/\__,_/\___/____/   \__,_/\___/_   \___/_/ /_/\____/\___/  

  ___  / /_   ____/ //_/ / /____  ____  / /____  _____                
 / _ \/ __/  / __  / _ \/ __/ _ \/ __ \/ __/ _ \/ ___/                
/  __/ /_   / /_/ /  __/ /_/  __/ / / / /_/  __(__  )                 
\___/\__/   \__,_/\___/\__/\___/_/ /_/\__/\___/____/                  
               __                          ___                        
  ____ ___  __/ /_____  __  _______   ____/ ( )__  ______  ___        
 / __ `/ / / / __/ __ \/ / / / ___/  / __  /|// / / / __ \/ _ \       
/ /_/ / /_/ / /_/ /_/ / /_/ / /     / /_/ /  / /_/ / / / /  __/       
\__,_/\__,_/\__/\____/\__,_/_/      \__,_/   \__,_/_/ /_/\___/        

  ____  ____ _(_)   _____                                             
 / __ \/ __ `/ / | / / _ \                                            
/ /_/ / /_/ / /| |/ /  __/                                            
\____/\__, /_/ |___/\___/                                             
     /____/                                                           
===============================================================================

Ondes de choc et détentes autour d'une ogive (Modèle 2D - plan)- Analyse complète
Créé le : 10 Décembre 2024
Auteur : MAFFIA Edoardo

Description :
Ce script calcule et analyse les propriétés thermodynamiques et aérodynamiques 
autour d'une ogive supersonique considérée ici comme plan. Les calculs incluent la pression statique, 
la pression totale, la température statique, la température totale, la vitesse locale,
les angles des ondes de choc et ondes de détendes.

Caractéristiques :
- Détermine les propriétés avant et après les ondes de choc ou les détentes de Prandtl-Meyer.
- Intègre des calculs pour la pression totale et la température totale à partir des relations isentropiques.
- Génère une illustration qui représente l'ogive et le tube avec les ondes de choc et de détente autour 
  de l'ogive, incluant les angles caractéristiques, les zones de transition (expansion Wave), et les lignes de Mach.

Applications :
- Comparaison et validation avec l'outil de simulation CFD Ansys.

Utilisation :
1. Définir les paramètres d'entrée : nombre de Mach initial, dimensions de l'ogive, 
   rapport des capacités thermiques (\(\gamma\)).
2. Exécuter le script pour obtenir les propriétés et visualiser les résultats graphiquement.

Remarque :
Les résultats fournis supposent un écoulement isentropique ou basé sur des relations classiques pour les ondes de choc.

ATTENTION:
CE MODELE N'EST UNIQUEMENT VALIDE SI L'ONDE DE CHOC EST ATTACHE (CHOC OBLIQUE)!
EN PLUS DE CA L'ETUDE CONSIDEREE EST PLAN ET NON AXYSYMETRIQUE

~~~~Pour des conditions réelles, des validations supplémentaires peuvent être nécessaires.~~~~

"""

#%%
# ===============================================================================
#     ____      __                 __           __  _           
#    /  _/___  / /__________  ____/ /_  _______/ /_(_)___  ____ 
#    / // __ \/ __/ ___/ __ \/ __  / / / / ___/ __/ / __ \/ __ \
#  _/ // / / / /_/ /  / /_/ / /_/ / /_/ / /__/ /_/ / /_/ / / / /
# /___/_/ /_/\__/_/   \____/\__,_/\__,_/\___/\__/_/\____/_/ /_/ 
# ===============================================================================
#Illustration pour différents Mach de la Thetha-Beta-Relation

from IPython.display import Image
from scipy.optimize import fsolve
from pylab import *
import matplotlib.pyplot as plt
import math
#%matplotlib inline

def func_theta_beta_mach(gamma, beta, mach):
    """
    Calcule l'angle de choc (theta) en fonction de :
    - l'exposant isentropique (gamma)
    - l'angle de déviation (delta)
    - le nombre de Mach (mach)
    """
    tan_theta = (2.0 * np.cos(beta) / np.sin(beta) *
                 (mach**2 * np.sin(beta)**2 - 1) /
                 (2.0 + mach**2 * (gamma + np.cos(2.0 * beta))))
    return np.arctan(tan_theta)

# --- Paramètres du problème ---
gamma = 1.4  # Exposant isentropique de l'air
beta_range = np.linspace(0.001, np.radians(90), num=80)  # beta de 0.001 rad à 90 deg (en radians)

# Liste de nombres de Mach pour lesquels on trace la relation
mach_numbers = np.array([1.1, 1.3, 1.5, 2.0, 2.5, 3.5, 4.8, 7.0, 15.0])

# Génération d'une palette de couleurs en dégradé de bleu
# Ici, on échantillonne la carte de couleurs 'Blues' de Matplotlib
colors = plt.cm.Blues(np.linspace(0.3, 0.9, len(mach_numbers)))
# Vous pouvez ajuster 0.3, 0.9 pour éclaircir ou foncer la palette

# --- Création de la figure et des axes ---
plt.figure(figsize=(8, 6))
plt.title("Relation de choc oblique", fontsize=16, fontweight='bold')

# Limites en degrés pour les angles de déviation (theta) et d'angle de choc (beta)
plt.xlim(0, 50)   # en degrés pour theta
plt.ylim(0, 90)   # en degrés pour beta

# --- Tracé des courbes pour chaque nombre de Mach ---
for i, M in enumerate(mach_numbers):
    theta_values = func_theta_beta_mach(gamma, beta_range, M)
    theta_degrees = np.degrees(theta_values)       # Conversion de l'angle de déviation en degrés
    beta_degrees = np.degrees(beta_range)          # Conversion de beta en degrés
    
    plt.plot(theta_degrees, 
             beta_degrees, 
             label=f"M = {M:.1f}", 
             linewidth=2,
             color=colors[i])  # Couleur tirée de la palette


# --- Personnalisation de l’affichage ---
plt.xlabel(r'Angle de déviation $\delta$ (deg)', fontsize=14)
plt.ylabel(r'Angle de choc $\beta$ (deg)', fontsize=14)
plt.grid(which='both', linestyle=':', linewidth=0.75, color='gray')
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

# Légende
plt.legend(title="Nombre de Mach", fontsize=12, title_fontsize=12)

# Ajustement de la mise en page
plt.tight_layout()

# --- Affichage ---
plt.show()



#%%
from IPython.display import Image
from scipy.optimize import fsolve
from pylab import *
import matplotlib.pyplot as plt
import math
#%matplotlib inline


def draw_cone_tube_horizontal(cone_length, cone_radius, tube_length, 
                              mach_angle=None, 
                              angle_mach_amont=None, 
                              angle_mach_aval=None):
    """
    Dessine une représentation 2D (en vue de coupe) d'un cône horizontal connecté à un tube.
    
    Paramètres obligatoires :
    -------------------------
    - cone_length : Longueur du cône
    - cone_radius : Rayon du cône à sa base
    - tube_length : Longueur du tube

    Paramètres optionnels :
    -----------------------
    - mach_angle       : Angle de l'onde de Mach en degrés (si spécifié, dessine les lignes d'onde de choc)
    - angle_mach_amont : Angle de la ligne de Mach "amont" (en degrés) partant du coin convexe
    - angle_mach_aval  : Angle de la ligne de Mach "aval" (en degrés) partant du coin convexe
    """
    
    # --------------------------------------------------------------------------------
    # 1) Définition des coordonnées du cône et du tube pour le tracé
    # --------------------------------------------------------------------------------
    # Le cône est défini par trois points (vu de profil) :
    #   - Sommet à x=0, y=0
    #   - Extrémité supérieure à x=cone_length, y=cone_radius
    #   - Extrémité inférieure à x=cone_length, y=-cone_radius
    cone_x = [0, cone_length, cone_length]
    cone_y = [0, cone_radius, -cone_radius]

    # Le tube est un rectangle (vu de profil) relié à la base du cône :
    #   - (cone_length, -cone_radius)
    #   - (cone_length + tube_length, -cone_radius)
    #   - (cone_length + tube_length, cone_radius)
    #   - (cone_length, cone_radius)
    tube_x = [cone_length, cone_length + tube_length, 
              cone_length + tube_length, cone_length]
    tube_y = [-cone_radius, -cone_radius, cone_radius, cone_radius]

    # --------------------------------------------------------------------------------
    # 2) Création de la figure et tracé des surfaces (cône + tube)
    # --------------------------------------------------------------------------------
    plt.figure(figsize=(8, 6))
    
    # On remplit la zone correspondant au cône
    plt.fill(cone_x, cone_y, color='gray', edgecolor='black', label='Cône + Tube')
    # On remplit la zone correspondant au tube
    plt.fill(tube_x, tube_y, color='gray', edgecolor='black')
    
    # --------------------------------------------------------------------------------
    # 3) Ajout de l’axe de symétrie (ligne pointillée)
    # --------------------------------------------------------------------------------
    symmetry_line_x = [-cone_length * 0.2, cone_length + tube_length]
    symmetry_line_y = [0, 0]
    plt.plot(symmetry_line_x, symmetry_line_y, color='black', 
             linestyle='--', label='Axe de symétrie')
    
    # --------------------------------------------------------------------------------
    # 4) Tracé de la prolongation du cône (lignes fines en pointillé)
    # --------------------------------------------------------------------------------
    cone_extension_x_top = [0, cone_length * 5]
    cone_extension_y_top = [0, cone_radius * 5]
    
    cone_extension_x_bottom = [0, cone_length * 5]
    cone_extension_y_bottom = [0, -cone_radius * 5]
    
    plt.plot(cone_extension_x_top, cone_extension_y_top, color='black', 
             linewidth=0.7, linestyle='--', label="Prolongement du cône")
    plt.plot(cone_extension_x_bottom, cone_extension_y_bottom, color='black', 
             linewidth=0.7, linestyle='--')
    
    # --------------------------------------------------------------------------------
    # 5) Calcul et affichage de l’angle du cône (delta) entre l’axe et la paroi
    # --------------------------------------------------------------------------------
    cone_angle_rad = np.arctan(cone_radius / cone_length)  # angle en radians
    cone_angle_deg = np.degrees(cone_angle_rad)            # angle en degrés
    
    # Arc de cercle pour représenter l’angle du cône
    angle_circle_radius = cone_length / 1.3
    theta = np.linspace(0, cone_angle_rad, 100)
    circle_x = angle_circle_radius * np.cos(theta)
    circle_y = angle_circle_radius * np.sin(theta)
    plt.plot(circle_x, circle_y, color='black')
    
    # Position du texte pour l'angle (delta)
    plt.text(circle_x[40] + 5, circle_y[40], 
             r'$\delta$ = ' + str(np.round(cone_angle_deg, 2)) + '°', 
             color='black', fontsize=10)

    # --------------------------------------------------------------------------------
    # 6) Tracé de l'onde de Mach (mach_angle) si spécifiée
    # --------------------------------------------------------------------------------
    if mach_angle is not None:
        mach_angle_rad = np.radians(mach_angle)
        
        # Lignes de l'onde de choc : on dessine deux droites symétriques par rapport à l'axe
        mach_line_x = [0, cone_length + tube_length]
        mach_line_y_top = [0, (cone_length + tube_length) * np.tan(mach_angle_rad)]
        mach_line_y_bottom = [0, -(cone_length + tube_length) * np.tan(mach_angle_rad)]
        
        plt.plot(mach_line_x, mach_line_y_top, color='red', linestyle='--', label='Onde de choc')
        plt.plot(mach_line_x, mach_line_y_bottom, color='red', linestyle='--')
        
        # Petit arc de cercle pour représenter l'angle (beta)
        circle_center_x = 0
        circle_center_y = 0
        circle_radius = cone_length / 5
        theta = np.linspace(0, mach_angle_rad, 100)
        circle_x = circle_center_x + circle_radius * np.cos(theta)
        circle_y = circle_center_y + circle_radius * np.sin(theta)
        plt.plot(circle_x, circle_y, color='red', linestyle='-', linewidth=2)
        
        # Position du texte pour l'angle de Mach (beta)
        plt.text(circle_x[70] + 5, circle_y[70], 
                 r'$\beta$ = ' + str(np.round(mach_angle, 2)) + '°', 
                 color='red', fontsize=12)

    # --------------------------------------------------------------------------------
    # 7) Tracé des lignes de Mach issues du coin convexe (angle_mach_amont / angle_mach_aval)
    # --------------------------------------------------------------------------------
    if (angle_mach_amont is not None) and (angle_mach_aval is not None):
        # Conversion degrés -> radians
        # On considère un "coin" à x=cone_length. 
        # Coin supérieur à y=cone_radius, coin inférieur à y=-cone_radius.
        angle_mach_aval_rad = np.radians(angle_mach_aval) 
        # Pour la ligne "amont", on ajoute l’angle du cône pour bien aligner avec la paroi
        angle_mach_amont_rad = np.radians(angle_mach_amont) + cone_angle_rad

        # Coin convexe supérieur (base du cône en haut)
        convex_corner_x_top = cone_length
        convex_corner_y_top = cone_radius

        # Droite aval (bleue) qui sort du coin supérieur
        extra_line1_x_top = [convex_corner_x_top, convex_corner_x_top + 600]
        extra_line1_y_top = [
            convex_corner_y_top,
            convex_corner_y_top + 600 * np.tan(angle_mach_aval_rad)
        ]
        plt.plot(extra_line1_x_top, extra_line1_y_top, color='blue', linestyle='--', 
                 label='Ligne de Mach aval')

        # Droite amont (verte) qui sort du coin supérieur
        extra_line2_x_top = [convex_corner_x_top, convex_corner_x_top + 600]
        extra_line2_y_top = [
            convex_corner_y_top,
            convex_corner_y_top + 600 * np.tan(angle_mach_amont_rad)
        ]
        plt.plot(extra_line2_x_top, extra_line2_y_top, color='green', linestyle='--', 
                 label='Ligne de Mach amont')

        # Coin convexe inférieur (base du cône en bas)
        convex_corner_x_bottom = cone_length
        convex_corner_y_bottom = -cone_radius

        # Droite aval (bleue) qui sort du coin inférieur
        extra_line1_x_bottom = [convex_corner_x_bottom, convex_corner_x_bottom + 600]
        extra_line1_y_bottom = [
            convex_corner_y_bottom,
            convex_corner_y_bottom + 600 * np.tan(-angle_mach_aval_rad)
        ]
        plt.plot(extra_line1_x_bottom, extra_line1_y_bottom, 
                 color='blue', linestyle='--')

        # Droite amont (verte) qui sort du coin inférieur
        extra_line2_x_bottom = [convex_corner_x_bottom, convex_corner_x_bottom + 600]
        extra_line2_y_bottom = [
            convex_corner_y_bottom,
            convex_corner_y_bottom + 600 * np.tan(-angle_mach_amont_rad)
        ]
        plt.plot(extra_line2_x_bottom, extra_line2_y_bottom, 
                 color='green', linestyle='--')

        # --------------------------------------------------------------------------------
        # 7a) Affichage de l’angle entre la ligne de Mach aval et l’axe horizontal (mu_2)
        # --------------------------------------------------------------------------------
        # On mesure l’angle entre la ligne aval (angle_mach_aval_rad) et l’horizontale (0 rad).
        angle_with_tube_rad = abs(angle_mach_aval_rad - 0)
        angle_with_tube_deg = np.degrees(angle_with_tube_rad)
        
        # Petit arc pour représenter l’angle mu_2
        angle_circle_radius = cone_length / 7
        theta = np.linspace(0, angle_mach_aval_rad, 100)
        circle_x = convex_corner_x_top + angle_circle_radius * np.cos(theta)
        circle_y = convex_corner_y_top + angle_circle_radius * np.sin(theta)
        plt.plot(circle_x, circle_y, color='blue')
        
        # Texte pour mu_2
        plt.text(circle_x[60] + 5, circle_y[60], 
                 r'$\mu_2$ = ' + str(np.round(angle_with_tube_deg, 2)) + '°', 
                 color='blue', fontsize=12)

        # --------------------------------------------------------------------------------
        # 7b) Affichage de l’angle entre la ligne de Mach amont et la prolongation du cône (mu_1)
        # --------------------------------------------------------------------------------
        # On définit la même logique mais pour l'angle amont (qui inclut l'angle du cône)
        prolongation_angle_rad = cone_angle_rad
        angle_between_rad = abs(angle_mach_amont_rad - prolongation_angle_rad)
        angle_between_deg = np.degrees(angle_between_rad)
        
        # Petit arc pour représenter l’angle mu_1
        angle_circle_radius = cone_length / 2
        theta = np.linspace(prolongation_angle_rad, angle_mach_amont_rad, 100)
        circle_x = convex_corner_x_top + angle_circle_radius * np.cos(theta)
        circle_y = convex_corner_y_top + angle_circle_radius * np.sin(theta)
        plt.plot(circle_x, circle_y, color='green')
        
        # Texte pour mu_1
        plt.text(circle_x[50] + 5, circle_y[50], 
                 r'$\mu_1$ = ' + str(np.round(angle_mach_amont, 2)) + '°', 
                 color='green', fontsize=12)
        
        angle_mach_aval_rad = np.radians(angle_mach_aval)
        angle_mach_amont_rad = np.radians(angle_mach_amont) + np.arctan(cone_radius / cone_length)

        angle_between_rad = abs(angle_mach_amont_rad - angle_mach_aval_rad)
        angle_between_deg = np.degrees(angle_between_rad)

        # Affichage de l'angle
        convex_corner_x_top = cone_length
        convex_corner_y_top = cone_radius
        angle_circle_radius = cone_length / 1.5
        theta = np.linspace(angle_mach_aval_rad, angle_mach_amont_rad, 100)
        circle_x = convex_corner_x_top + angle_circle_radius * np.cos(theta)
        circle_y = convex_corner_y_top + angle_circle_radius * np.sin(theta)
        plt.plot(circle_x, circle_y, color='purple', linestyle='-')
        plt.text(circle_x[99], circle_y[99]+10, r'$\gamma$=' + str(np.round(angle_between_deg, 2)) + '°', color='purple', fontsize=12)
    
    # --------------------------------------------------------------------------------
    # 8) Paramètres d'affichage final
    # --------------------------------------------------------------------------------
    plt.axis('equal')  # Même échelle sur X et Y
    plt.xlim(-cone_length * 0.2, cone_length + tube_length)  
    plt.ylim(-cone_radius * 2, cone_radius * 2)  
    plt.title("Onde de Mach & Détente de Prandtl-Meyer")
    plt.xlabel("Longueur (mm)")
    plt.ylabel("Hauteur (mm)")
    plt.grid(True)
    plt.legend(fontsize=8)
    plt.show()

# ------------------------------------------------------------------------------
# 1) Définition de constantes et de paramètres
# ------------------------------------------------------------------------------
gam = 1.4              # Exposant isentropique (ratio des chaleurs spécifiques)

h = 405                # Longueur du cône (mm)
r = 98/2               # Rayon du cône (mm) -> 98 mm de diamètre, donc 49 mm de rayon
l = 700                # Longueur du tube (mm)

# ------------------------------------------------------------------------------
# 4) Tracé de la géométrie (cône + tube) via la fonction draw_cone_tube_horizontal
# ------------------------------------------------------------------------------
draw_cone_tube_horizontal(cone_length=h, cone_radius=r, tube_length=l)
delta = np.arcsin(r/h) * 180.0 / np.pi

#%%
# Liste de nombres de Mach pour lesquels on trace la relation
mach_numbers_test = np.array([1.3, 1.31, 1.32, 1.33, 1.34, 1.35])

beta_range = np.linspace(0.001, np.radians(90), num=80)  # beta de 0.001 rad à 90 deg (en radians)

# Génération d'une palette de couleurs en dégradé de bleu
# Ici, on échantillonne la carte de couleurs 'Blues' de Matplotlib
colors = plt.cm.Blues(np.linspace(0.3, 0.9, len(mach_numbers_test)))
# Vous pouvez ajuster 0.3, 0.9 pour éclaircir ou foncer la palette

# --- Création de la figure et des axes ---
plt.figure(figsize=(8, 6))
plt.title("Relation de choc oblique", fontsize=16, fontweight='bold')

# Limites en degrés pour les angles de déviation (theta) et d'angle de choc (beta)
plt.xlim(0, 10)   # en degrés pour theta
plt.ylim(40, 90)   # en degrés pour beta

# --- Tracé des courbes pour chaque nombre de Mach ---
for i, M_test in enumerate(mach_numbers_test):
    theta_values = func_theta_beta_mach(gam, beta_range, M_test)
    theta_degrees = np.degrees(theta_values)       # Conversion de l'angle de déviation en degrés
    beta_degrees = np.degrees(beta_range)          # Conversion de beta en degrés
    
    plt.plot(theta_degrees, 
             beta_degrees, 
             label=f"M = {M_test:.2f}", 
             linewidth=2,
             color=colors[i])  # Couleur tirée de la palette

plt.axvline(x=delta, color='gray', linestyle='--', alpha=0.5)

# --- Personnalisation de l’affichage ---
plt.xlabel(r'Angle de déviation $\delta$ (deg)', fontsize=14)
plt.ylabel(r'Angle de choc $\beta$ (deg)', fontsize=14)
plt.grid(which='both', linestyle=':', linewidth=0.75, color='gray')
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

# Légende
plt.legend(title="Nombre de Mach", fontsize=12, title_fontsize=12)

# Ajustement de la mise en page
plt.tight_layout()

# --- Affichage ---
plt.show()

#%%

# ------------------------------------------------------------------------------
# 1) Définition des paramètres
# ------------------------------------------------------------------------------
Ma1 = 1.33         # On prend le Mach idéal
P1 = 90178.8  # Pression en (Pa) -> à la hauteur 1000m (vitesse maximale)
T1 = 13.5 + 273.15  # Température en (K) -> à la hauteur 1000m (vitesse maximale)

# ------------------------------------------------------------------------------
# 3) Calcul des grandeurs stagnation (Po1, To1) en conditions isentropiques
#    relations isentropiques en écoulement compressible d’un gaz parfait
# ------------------------------------------------------------------------------
Po1 = P1 * (1 + (gam - 1) / 2 * Ma1**2)**(gam / (gam - 1))
To1 = T1 * (1 + (gam - 1) / 2 * Ma1**2)



# ------------------------------------------------------------------------------
# 5) affichage des résultats

print(f"δ = {delta:.2f}°")
print("Condition en Amont:")
print("Ma1 = ",Ma1)
print(f"P1  = {P1} Pa")
print(f"T1  = {T1} K")
print(f"Po1 = {Po1} Pa")
print(f"To1 = {To1} K")
print()

#%%

# ----------------------------------------------------------------------
# 1) Définition de la fonction qui calcule l'angle de déviation (theta)
#    en fonction de gamma (gam), beta et du Mach M
# ----------------------------------------------------------------------
def func_Theta_Beta_Mach(gam, beta, M):
    """
    Calcule l'angle de déviation (theta) pour un choc oblique
    en fonction de:
      - gam : coefficient de chaleur spécifique (1.4 pour l'air)
      - beta : angle de choc (shock angle) en radians
      - M : nombre de Mach en amont du choc.

    Formule basée sur la relation de l'onde de choc oblique :
        tan(theta) = 2 * cot(beta) * ((M^2 sin^2(beta) - 1) / [2 + M^2 (gam + cos(2 beta))])
    
    Retourne : theta (en radians)
    """
    # Expression du tan(theta)
    numerator = 2.0 * cos(beta) / sin(beta) * (M**2 * sin(beta)**2 - 1)
    denominator = (2.0 + M**2 * (gam + cos(2.0 * beta)))
    tan_theta = numerator / denominator
    
    # Conversion en angle (radians)
    theta = arctan(tan_theta)
    return theta


# ----------------------------------------------------------------------
# 2) Fonction pour trouver les angles de choc (beta) qui correspondent
#    à un angle de déviation (delta_target) donné.
#    Utilise la foction fsolve de la bibliothèque SciPy
#    (méthode basée sur des techniques dérivées de Newton-Raphson)
# ----------------------------------------------------------------------
def find_intersections(gam, M, delta_target):
    """
    Recherche numérique (via fsolve) des solutions "faible" et "forte"
    pour beta (en degrés), telles que theta(beta) = delta_target.

    - gam          : rapport des chaleurs spécifiques
    - M            : nombre de Mach
    - delta_target : angle de déviation visé (en degrés)

    Retourne : (beta_weak, beta_strong) en degrés
    """

    # Définition de l'équation à résoudre : f(beta) = theta_deg - delta_target = 0
    def equation(beta):
        # beta est fourni en degrés, conversion en radians pour le calcul
        beta_rad = beta * np.pi / 180.0
        theta = func_Theta_Beta_Mach(gam, beta_rad, M)
        theta_deg = theta * 180.0 / np.pi
        return theta_deg - delta_target
    
    # Recherche d'une solution "faible" (weak shock), initial guess à 20°
    beta1 = fsolve(equation, 20)[0]
    
    # Recherche d'une solution "forte" (strong shock), initial guess à 80°
    beta2 = fsolve(equation, 80)[0]
    
    return beta1, beta2


# ----------------------------------------------------------------------
# 3) Paramètres généraux pour la démonstration et le tracé
# ----------------------------------------------------------------------

beta_range = np.linspace(0.001, 90, num=80) * pi / 180.0 
# beta_range : angles de choc de ~0.001° à 90°, convertis en radians.


# Calcul de l'angle du cône (delta) en degrés : arcsin(r/h)
delta = np.arcsin(r / h) * 180.0 / np.pi

# ----------------------------------------------------------------------
# 4) Tracé de la relation choc oblique (theta-beta) pour Mach = M
# ----------------------------------------------------------------------
plt.figure(figsize=(8, 6))
plt.xlim(0, 10)   # Limites sur l'axe X : angle de déviation (delta) en degrés
plt.ylim(40, 90)   # Limites sur l'axe Y : angle de choc (beta) en degrés


# On calcule theta (en radians) pour la plage beta_range (en radians),
# puis on convertit theta en degrés.
theta_deg = func_Theta_Beta_Mach(gam, beta_range, Ma1) * 180.0 / pi
    
# On trace beta (en degrés) en fonction de theta (en degrés).
plt.plot(theta_deg, beta_range * 180.0 / pi, label=f"M={Ma1:.1f}", linewidth=2)

# Ajouts de titre, légendes, grilles, etc.
plt.title("Oblique Shock Relation", fontsize=14)
plt.xlabel(r'Angle de déviation $\delta$ (deg)', fontsize=14)
plt.ylabel(r'Angle de Mach $\beta$ (deg)', fontsize=14)
plt.grid(which='both', linestyle='-', linewidth=0.75, color='0.85')
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

# ----------------------------------------------------------------------
# 5) Recherche des points d'intersection (weak & strong shocks)
#    pour M = 2 et l'angle de cône delta
# ----------------------------------------------------------------------
beta_weak, beta_strong = find_intersections(gam, Ma1, delta)

print(f"\nPoints d'intersection pour M={Ma1} à δ = {delta:.2f}° :")
print(f"  -> Solution faible (weak shock)  : β = {beta_weak:.2f}°")
print(f"  -> Solution forte (strong shock) : β = {beta_strong:.2f}°")

# ----------------------------------------------------------------------
# 6) Ajout des points d'intersection sur le graphique
# ----------------------------------------------------------------------
# On place un marqueur pour la solution faible (rouge) et la solution forte (bleue).
plt.plot(delta, beta_weak, 'ro', label='Solution faible', markersize=8)
plt.plot(delta, beta_strong, 'bo', label='Solution forte', markersize=8)

# On trace des lignes pointillées pour repérer ces points sur le graphique
plt.axvline(x=delta, color='gray', linestyle='--', alpha=0.5)
plt.axhline(y=beta_weak, color='red', linestyle=':', alpha=0.5)

plt.legend(title="Mach Number")
plt.tight_layout()
plt.show()


# ----------------------------------------------------------------------
# 7) Tracé de la géométrie (cône + tube) en utilisant l'angle "faible"
#    en tant qu'onde de Mach (mach_angle)
# ----------------------------------------------------------------------
draw_cone_tube_horizontal(
    cone_length=h, 
    cone_radius=r, 
    tube_length=l, 
    mach_angle=beta_weak
)
print(f"\n")
#%%

# ----------------------------------------------------------------------
# 1) Calcul du Mach normal en amont du choc (Mn1)
#    beta_weak est en degrés, on le convertit donc en radians (π/180).
# ----------------------------------------------------------------------
Mn1 = Ma1 * np.sin(beta_weak * np.pi / 180.0)
print("Vitesse Normale:")
print(f"M_n1 = {Mn1:.3f}")

# ----------------------------------------------------------------------
# 2) Calcul du Mach normal en aval du choc (Mn2)
#    Relation du choc normal : 
#    Mn2 = sqrt( ((γ - 1) * Mn1² + 2) / (2γ * Mn1² - (γ - 1)) )
# ----------------------------------------------------------------------
Mn2 = sqrt(((gam - 1.0) * Mn1**2 + 2.0) / (2.0 * gam * Mn1**2 - (gam - 1.0)))
print(f"M_n2 = {Mn2:.3f}")

# ----------------------------------------------------------------------
# 1) Calcul du Mach en aval du choc (Ma2)
#    beta_weak est en degrés, on le convertit donc en radians (π/180).
# ----------------------------------------------------------------------
Ma2 = Mn2 / sin((beta_weak - delta)* np.pi/180)
print(f"Ma2  = {Ma2:.3f}")

print(f"\n")
#%%
# --------------------------------------------------------------
# Relations de choc normal (Adiabatique et réversible) :
# --------------------------------------------------------------
# 1) Pression statique en aval (P2)
P2 = ( (2.0 * gam * Mn1**2 - (gam - 1.0)) / (gam + 1.0) ) * P1


# 2) Température statique en aval (T2)
T2 = ( ((2.0 * gam * Mn1**2 - (gam - 1.0)) * ((gam - 1.0) * Mn1**2 + 2.0)) 
       / ((gam + 1.0)**2 * Mn1**2) ) * T1

# 3) Pression totale en aval (Po2)
Po2 = ( ((gam + 1.0)* Mn1**2) / ((gam - 1.0)* Mn1**2 + 2.0) )**( gam / (gam - 1.0) ) \
      * ( ((gam + 1.0) / (2.0 * gam * Mn1**2 - (gam - 1.0))) )**( 1.0 / (gam - 1.0) ) \
      * Po1

# 4) Température totale en aval (To2)
#    Dans un choc normal, la température totale (adiabatique) est conservée.
To2 = To1

# --------------------------------------------------------------
# Affichage des résultats
# --------------------------------------------------------------
print("Choc normal :")
print(f"P2   = {P2:.2f} Pa")
print(f"T2   = {T2:.2f} K")
print(f"Po2  = {Po2:.2f} Pa")
print(f"To2  = {To2:.2f} K")

print(f"\n")

#%%

import math
from scipy.optimize import fsolve


# ------------------------------------------------------------------------------
# Fonction de Prandtl-Meyer
# ------------------------------------------------------------------------------
def prandl_meyer(gam, M):
    """
    Calcule la fonction de Prandtl-Meyer ν(M) (en radians) pour un écoulement supersonique.
    gam : rapport des chaleurs spécifiques (ex. 1.4 pour l'air)
    M   : nombre de Mach

    Formule :
        ν(M) = sqrt((gam + 1)/(gam - 1)) * arctan( sqrt(((gam - 1)/(gam + 1)) * (M^2 - 1)) )
               - arctan( sqrt(M^2 - 1) )
    Retourne : ν(M) en radians
    """
    nu = np.sqrt((gam + 1)/(gam - 1)) * np.arctan(np.sqrt(((gam - 1)/(gam + 1)) * (M**2 - 1))) - np.arctan(np.sqrt(M**2 - 1))
    return nu


# ------------------------------------------------------------------------------
# 1) Calcul de l'angle ν(Ma2) en degrés
# ------------------------------------------------------------------------------
nu2 = prandl_meyer(gam, Ma2)  # en radians
nu2_deg = np.degrees(nu2)     # conversion en degrés

print("Onde de détente")

print("ν(Ma2)  = ", nu2_deg, "°")

# ------------------------------------------------------------------------------
# 2) Calcul de nu3 = delta + ν(Ma2)
#    Ici, on ajoute l'angle delta (déviation supplémentaire) à l'angle ν(Ma2).
#    ATTENTION : delta doit être en degrés ; on additionne donc en degrés.
# ------------------------------------------------------------------------------
nu3 = delta + nu2_deg
print("ν(Ma3)  = ", nu3, "°")

# ------------------------------------------------------------------------------
# 3) Définition de la fonction f(Ma3) = ν(Ma3) - ν3
#    On veut f(Ma3) = 0  =>  ν(Ma3) = ν3
#    => solve for Ma3
# ------------------------------------------------------------------------------
def f(Ma3, gam, nu3_deg):
    """
    Évalue la différence : ν(Ma3) - nu3.
    Ma3       : nombre de Mach inconnu
    gam       : gamma
    nu3_deg   : valeur cible de la fonction Prandtl-Meyer (en degrés)

    Retourne f(Ma3) = ν(Ma3) (en degrés) - nu3_deg
    """
    # ν(Ma3) en radians
    nu_Ma3 = prandl_meyer(gam, Ma3)
    # conversion en degrés
    nu_Ma3_deg = nu_Ma3 * 180/np.pi

    return nu_Ma3_deg - nu3_deg

# ------------------------------------------------------------------------------
# 4) Résolution numérique via fsolve pour trouver Ma3
#    - Initial guess : 1.2 (doit être >1 pour écoulement supersonique)
#    Utilise la foction fsolve de la bibliothèque SciPy
#    (méthode basée sur des techniques dérivées de Newton-Raphson)
# ------------------------------------------------------------------------------
Ma3_guess = 1.1
Ma3_solution = fsolve(f, Ma3_guess, args=(gam, nu3))[0]  # on prend le premier élément du tableau fsolve

print("Ma3     = ", Ma3_solution)

# ------------------------------------------------------------------------------
# 5) Angles de Mach (µ1, µ2) = arcsin(1 / M) en radians, puis conversion en degrés
# ------------------------------------------------------------------------------
mu1_rad = np.arcsin(1.0 / Ma2)
mu2_rad = np.arcsin(1.0 / Ma3_solution)

mu1_deg = mu1_rad * 180/np.pi
mu2_deg = mu2_rad * 180/np.pi

print("µ1      = ", mu1_deg, "°")
print("µ2      = ", mu2_deg, "°")

draw_cone_tube_horizontal(
    cone_length=h, 
    cone_radius=r, 
    tube_length=l, 
    mach_angle=beta_weak,
    angle_mach_amont=mu1_deg-delta,
    angle_mach_aval=mu2_deg
)
#ATTENTION DE BIEN SOUSTRAIRE PAR L'ANGLE DE DEVIATION DELTA POUR REPRESENTER L'ANGLE EN AMONT!!!
print(f"\n")

#%%
# Pression et température totales conservées (onde de détente isentropique)
Po3 = Po2  
To3 = To2  

# Calcul de P3 et T3 via les formules d'écoulement isentropiques
T3 = T2 * (1 + (gam - 1)/2 * Ma2**2) / (1 + (gam - 1)/2 * Ma3_solution**2)
P3 = P2 * ( (T3 / T2) ** (gam / (gam - 1)) )

print("P3  =", P3, "Pa")
print("T3  =", T3, "K")
print("Po3 =", Po3, "Pa")
print("To3 =", To3, "K")

