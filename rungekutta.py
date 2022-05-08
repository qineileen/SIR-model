import sympy

class RungeKutta():
    """
    this is for rungekutta method
    """
    def __init__(self, f, h, x0, y0, a, b):
        """
        f is funcion
        h is step size
        x0, y0 are initial value
        [a , b] is search range
        n is loop number
        """
        self.f = f
        self.h = h
        self.x0 = x0
        self.y0 = y0
        self.n = int((b-a)/h)

    def solve_function1(self):
        """
        dy/dx= 1/(x**2)-y/x
        """

        Y = []
        X = []
        Y.append(self.y0)
        X.append(self.x0)
        
        for i in range(self.n):
            k1 = self.f.subs([(x, self.x0), (y, Y[i])])
            k2 = self.f.subs([(x, self.x0+self.h/2), (y, Y[i]+self.h*k1/2)])
            k3 = self.f.subs([(x, self.x0+self.h/2), (y, Y[i]+self.h*k2/2)])
            k4 = self.f.subs([(x, self.x0+self.h), (y, Y[i]+self.h*k3)])
            y_rk4 = Y[i]+self.h*(k1+2*k2+2*k3+k4)/6
            self.x0 = self.x0+self.h
            Y.append(y_rk4)
            X.append(self.x0)
        return X, Y


class Euler:
    """
    f is funcion
        h is step size
        x0, y0 are initial value
        [a , b] is search range
        n is loop number
    """
    def __init__(self, f, h, x0, y0, a, b): 
        self.f = f
        self.h = h
        self.x0 = x0
        self.y0 = y0
        self.n = int((b-a)/h)

    def solve_function2(self):
        X=[]
        Y=[]
        X.append(self.x0)
        Y.append(self.y0)
        for i in range(self.n):
            y_euler=Y[i]+self.h*self.f.subs([(x,self.x0),(y,Y[i])])
            self.x0=self.x0+self.h
            Y.append(y_euler)
            X.append(self.x0)
        return X, Y


def solve_RungeKutta(x, y):
    f = 1/(x**2)-y/x
    t1 = RungeKutta(f, 0.005, 1, 1, 1, 2)
    x1, y1 = t1.solve_function1()
    for i in range(len(x1)):
        print("xi:", "%.3f" % x1[i], "     y_rk4:", "%.7f" % y1[i])

def solve_Euler(x, y):
    f=y-2*x/y
    t2 = Euler(f,0.1,0,1,0,1)
    x2, y2 = t2.solve_function2()
    for i in range(len(x2)):
        print("xi:", "%.3f" % x2[i], "     y_euler:", "%.7f" % y2[i])




class ProEuler:
     def __init__(self, f, h, x0, y0, a, b):
         self.h=h
         self.x0=x0
         self.y0=y0
         self.f=f
         self.n=int((b-a)/h)

     def solve_function3(self):
         X=[]
         Y=[]
         X.append(self.x0)
         Y.append(self.y0)
         for i in range(self.n):
             y_euler1=Y[i]+self.h*self.f.subs([(x,self.x0),(y,Y[i])])
             self.x0=self.x0+self.h
             y_euler2=Y[i]+self.h*self.f.subs([(x,self.x0),(y,y_euler1)])
             y_pro_euler=(y_euler1+y_euler2)/2
             Y.append(y_pro_euler)
             X.append(self.x0)
         return X, Y

def solve_ProEuler(x,y):
    f=y-2*x/y
    t3=ProEuler(f,0.1,0,1,0,1)
    x3, y3 = t3.solve_function3()
    for i in range(len(x3)):
        print("xi:", "%.3f" % x3[i], "     y_proeuler:", "%.7f" % y3[i])




if __name__ == "__main__":
    x, y = sympy.symbols("x  y")
    #solve_Euler(x, y)
    #solve_RungeKutta(x, y)
    solve_ProEuler(x,y)
    
