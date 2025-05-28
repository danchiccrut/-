import urllib.request
import math
import os
import re
import scipy.special as scp
import matplotlib.pyplot as plt


class Rcs:
    def __init__(self, r, fmin, fmax):
        self.r = r
        self.fmin = fmin
        self.fmax = fmax
        self.mass_f = []
        self.mass_rcs = []

    def calc(self):
        def hn(n, x):
            return complex(scp.spherical_jn(n, x), scp.spherical_yn(n, x))

        f = self.fmin
        print(self.fmin, self.fmax)
        while f <= self.fmax:
            l = 3 * 10**8 / f
            k = 2 * math.pi / l
            s = 0
            for n in range(1, 20):
                an = scp.spherical_jn(n, k * self.r) / hn(n, k * self.r)
                bn = ((k * self.r * scp.spherical_jn(n - 1, k * self.r) - n * scp.spherical_jn(n, k * self.r)) /
                      (k * self.r * hn(n - 1, k * self.r) - n * hn(n, k * self.r)))
                s += ((-1) ** n) * (n + 0.5) * (bn - an)
            self.mass_f.append(f)
            self.mass_rcs.append((l**2) * ((abs(s))**2) / math.pi)
            f += 40_000_000

        plt.plot(self.mass_f, self.mass_rcs)
        plt.xlabel("f, Hz")
        plt.ylabel("RCS, m^2")
        plt.grid(True)
        plt.show()


class Output:
    def __init__(self, mass_f, mass_rcs):
        self.mass_f = mass_f
        self.mass_rcs = mass_rcs

    def output(self):
        results_dir = "results"
        if not os.path.isdir(results_dir):
            os.mkdir(results_dir)

        file_path = os.path.join(results_dir, "res2.xml")

        with open(file_path, 'w', encoding='utf-8') as f:
            f.write('<data>\n\t')
            f.write('<frequencydata>\n\t\t')
            for i in self.mass_f:
                f.write('<f>{:e}</f>\n\t\t'.format(i))
            f.write('</frequencydata>\n\t')
            f.write('<lambdadata>\n\t\t')
            for i in self.mass_f:
                l = 3 * 10**8 / i
                f.write('<lambda>{:e}</lambda>\n\t\t'.format(l))
            f.write('</lambdadata>\n\t')
            f.write('<rcsdata>\n\t\t')
            for i in self.mass_rcs:
                f.write('<rcs>{:e}</rcs>\n\t\t'.format(i))
            f.write('</rcsdata>\n\t')
            f.write('</data>\n')


# Загрузка входного файла
url = 'https://jenyay.net/uploads/Student/Modelling/task_rcs_01.csv'
urllib.request.urlretrieve(url, "input.csv")

# Извлекаем параметры из строки, начинающейся на "1"
with open("input.csv", "r", encoding="utf-8") as F:
    for line in F:
        if line.startswith("1"):
            numbers = re.findall(r'\d+(?:\.\d+)?(?:[eE][-+]?\d+)?', line)
            numbers = [float(n) for n in numbers]
            task_num, fmin, fmax, D = numbers[:4]
            r = D / 2
            break

# Запуск расчётов и сохранения
sph = Rcs(r, fmin, fmax)
sph.calc()

o = Output(sph.mass_f, sph.mass_rcs)
o.output()

