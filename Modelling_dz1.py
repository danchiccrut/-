import os
import csv
import matplotlib.pyplot as plt
from math import*

#Определение переменных и функции
A=10
def y(x):
    return A+x**2-cos(2*pi*x)
X=[]
Y=[]

#создание директории 

if not os.path.isdir("results"): 
    os.makedirs("results") 

#сoхраняем в csv

c=-5.12 

while (c<5.13): 
    X.append(c) 
    Y.append(y(c)) 
    c+=0.1 

with open('results/result.csv', 'w', newline='', encoding='utf-8') as f:
    writer = csv.writer(f)
    
    #запишем название столбцов
    writer.writerow(['X', 'Y'])
    print(f"{'X'}{'Y':>7}")
    
    #запишем пары значений по строкам
    for x_val, y_val in zip(X, Y):
        writer.writerow([f"{x_val:.3f}", f"{y_val:.3f}"])
        print(f"{x_val:.3f}", f"{y_val:.3f}")




f.close() 
plt.plot(X,Y) 
plt.show()
