#Calculates the pillowcase homology for Montesinos knots
import numpy as np
import time
from matplotlib.ticker import FuncFormatter
import snappy #Used to compute signatures

#Makes the axes label on sage graphs in terms of pi
def f(x,pos):
    if x == 0:
        return r'0'
    elif x==1:
        return r'$\pi$'
    elif x==2:
        return r'$2\pi$'
    elif x==1/2:
        return r'$\frac{\pi}{2}$'
    elif x==3/2:
        return r'$\frac{3\pi}{2}$'
    else:
        return r'\?'
formatter = FuncFormatter(f)


class Tangle:
    def __init__(self, paths):
        self.paths = paths

    def breakall(self, height):
        newpaths = []
        for p in self.paths:
            if (p.start[1]<height and p.end[1]>height) or (p.start[1]>height and p.end[1]<height):
                slope = (p.start[1]-p.end[1])/(p.start[0]-p.end[0])
                x = p.start[0]+(height - p.start[1])/slope
                p1 = Path(p.start,(x,height),p.prev)
                p2 = Path((x,height),p.end,p1,p.next)
                p1.next = p2
                if p.prev != None:
                    p.prev.next=p1
                if p.next != None:
                    p.next.prev=p2
                newpaths.append(p1)
                newpaths.append(p2)
            else:
                newpaths.append(p)
        self.paths = newpaths

    def getstart(self, point):
        for p in self.paths:
            if p.start == point or (p.start[0],p.start[1]-2)==point:
                return p

    def getend(self, point):
        for p in self.paths:
            if p.end == point or (p.end[0],p.end[1]-2)==point:
                return p

    def split(self):
        #Split into multiple components (accounting for perturbation)
        #ps = self.paths.copy()
        #TODO: Does not work if earring has been done first
        ps = [x for x in self.paths]
        curr = self.getstart((0,0))
        connpath = [curr]
        ps.remove(curr)
        while curr.next != None:
            #print(curr)
            curr = curr.next
            if curr.next == connpath[-1]:
                curr.reverse()
            connpath.append(curr)
            ps.remove(curr)
        #print(len(connpath))
        #components = [Tangle(deepcopy(connpath))]
        components = [Tangle(copy(connpath))]

        #used = connpath.copy()

        while len(ps)>0:
            curr = ps[0]
            currpath = [curr]
            ps.remove(curr)
            while curr.next in ps:
                curr = curr.next
                if curr.next == currpath[-1]:
                    curr.reverse()
                currpath.append(curr)
                ps.remove(curr)
            #components.append(Tangle(deepcopy(currpath)))
            components.append(Tangle(copy(currpath)))

        return(components)

    def invert(self):
        #Reflects image across the line x=y. Corresponds to change of parametrization from reflecting tangle across line x=y
        pass

    def normx(self,x):
        if x<0:
            return -x
        if x>1:
            return 2-x
        return x

    def show(self,N = True):
        show(self.getimg(),aspect_ratio=1,figsize=12, tick_formatter=[formatter,formatter], ticks=[1/2,1/2], axes_labels=[r'$\gamma$',r'$\theta$'])

    def getcompimg(self,width = 300, N=True, buff = 10):
        comps = self.split()
        n = len(comps)
        im = np.ones((width*2+buff*2,width*n+buff*(n+1),4),np.uint8)*0
        for i, c in enumerate(comps):
            im[buff:buff+2*width, (i+1)*buff+i*width:(i+1)*buff+(i+1)*width] = c.getimg(N,width)
        return im


    def getimg(self,flaglines=True,colup='blue',coldown='blue'):
        im = Graphics()
        im.axes_range(0,1,0,2)
        im += line([(0,1),(1,1)], rgbcolor=(0,0,0),thickness = 1)
        n = max([l.flag for l in self.paths])+1
        #nf = 5 #Avoids blue and yellow
        colors = rainbow(n)
        for p in self.paths:
            if p.flag:
                if flaglines:
                    im += line([p.start,p.end],color=colup,thickness = 2)
                else:
                    im += circle(p.start,0.1,rgbcolor = colors[p.flag],fill = True)
                    im += circle(p.end,0.1,rgbcolor = colors[n],fill = True)
                #nf += 1
            else:
                if p.end[1]>=p.start[1]:
                    im += line([p.start,p.end],rgbcolor=colup,thickness=2)
                else:
                    im += line([p.start,p.end],rgbcolor=coldown,thickness=2)

        return im


    def getcorner(self, paths, last = False):
        for p in paths:
            if p.prev == None and not last:
                return p
            if p.next == None and last:
                return p

    def earring(self, d=1/10,smart=True):
        #NOTE:Only works correctly on the arc component
        if smart:
            if abs(self.paths[0].slope) > 1:
                f = self.getcorner(self.paths)
                if f.slope>0:
                    self.subdivide(f,(d/f.slope,d))
                else:
                    self.subdivide(f,(-d/f.slope,2-d)) 
                b = self.getcorner(self.paths,True)
                if b.end[0] == 1 and b.slope>0:
                    self.subdivide(b,(1-d/b.slope,b.end[1]-d))
                elif b.end[0] == 1 and b.slope<0:
                    self.subdivide(b,(1+d/b.slope,b.end[1]+d))
                elif b.end[0] == 0 and b.slope>0:
                    self.subdivide(b,(d/b.slope,b.end[1]+d))
                else:
                    self.subdivide(b,(-d/b.slope,b.end[1]-d))
            else:
                f = self.getcorner(self.paths)
                self.subdivide(f,(d,f.start[1]+f.slope*d))
                b = self.getcorner(self.paths,True)
                if b.end[0] == 0:
                    self.subdivide(b,(d,b.end[1]+b.slope*d))
                else:
                    self.subdivide(b,(1-d,b.end[1]-b.slope*d))
        d/=5/2
        pathsf = self.paths.copy()
        pathsb = deepcopy(self.paths)
        conn = []
        for p in pathsb:
            p.reverse()
        fs = self.getcorner(paths = pathsf)
        fe = self.getcorner(pathsf, True)
        bs = self.getcorner(paths = pathsb)
        be = self.getcorner(pathsb, True)
        dp = [fs,be] #Paths that might contain the distinguished intersection

        if fs.start == (0,0):
            fs.start = (0,d)
            be.end = (d,0)
            p = Path((d,2),(0,2-d),be,fs)
        elif fs.start == (0,2):
            fs.start = (d,2)
            be.end = (0,2-d)
            p = Path((0,d),(d,0),be,fs)
        fs.prev = p
        be.next = p
        conn.append(p)
        dp.append(p)

        if fe.end == (1,0):
            fe.end = (1-d,0)
            bs.start = (1,d)
            p = Path((1-d,2),(1,2-d),fe,bs)
            bs.prev = p
            fe.next = p
            conn.append(p)
        elif fe.end == (1,2):
            fe.end = (1,2-d)
            bs.start = (1-d,2)
            p = Path((1,d),(1-d,0),fe,bs)
            bs.prev = p
            fe.next = p
            conn.append(p)
        elif fe.end == (0,1):
            fe.end = (0,1+d)
            bs.start = (0,1-d)
            bs.prev = fe
            fe.next = bs
        elif fe.end == (1,1):
            fe.end = (1,1-d)
            bs.start = (1,1+d)
            bs.prev = fe
            fe.next = bs

        self.paths = pathsf + pathsb + conn
        for p in self.paths:
            p.reslope()
        return dp

    def subdivide(self,path,coord):
        #subdivide a path at a particular coordinate (usefull for intersections)
        if not equiv(coord,path.start) and not equiv(coord,path.end):
            if path.start[1]==0 and path.end[1]==0 and coord[1] == 2:
                p1 = Path(path.start,(coord[0],0),path.prev,flag = path.flag)
                p2 = Path(coord,(path.end[0],2),p1,path.next,path.flag)
            elif path.start[1]==2 and path.end[1]==2 and coord[1] == 0:
                p1 = Path(path.start,(coord[0],2),path.prev,flag = path.flag)
                p2 = Path(coord,(path.end[0],0),p1,path.next,path.flag)
            else:
                p1 = Path(path.start,coord,path.prev,flag = path.flag)
                p2 = Path(coord,path.end,p1,path.next,path.flag)
            p1.next = p2
            if path.prev != None:
                path.prev.next = p1
            if path.next != None:
                path.next.prev = p2
            i = self.paths.index(path)
            self.paths.remove(path)
            self.paths.insert(i, p1)
            self.paths.insert(i+1, p2)
            return [p1,p2]
        return None

    def showlift(self):
        #Show the lift of the character variety.
        p = self.getcorner(self.paths)
        paths = [p.start,p.end]
        while p.next != None:
            p=p.next
            paths.append(p.start)
            paths.append(p.end)
        i = IntersectionPoint((0,0),None,None)
        b = Bigon(i,i,paths)
        show(line(b.lift))

    #TODO super getstart, getend, getc0, getcpi

class RatTangle(Tangle):
    def __init__(self,slope):
        self.slope = slope
        self.paths = []
        self.gen()

    def gen(self):
        self.paths = [None]#Hacky
        direction = True # self.slope>0 #True if orientation is to the right
        if self.slope>0:
            curr = (0,0)
        else:
            curr = (0,2)
        endpoints = [(0,0),(0,1),(1,0),(1,1),(0,2),(1,2)]
        endpoints.remove(curr)
        while curr not in endpoints:
            if direction:
                (x,y) = (1,curr[1] + self.slope*(1-curr[0]))
                if y <= 2 and y>=0:
                    self.paths.append(Path(curr,(x,y),self.paths[-1]))
                    curr = (x,2-y)
                    direction = not direction
                elif self.slope > 0:
                    (x,y) = (curr[0]+(2-curr[1])/self.slope,2)
                    self.paths.append(Path(curr,(x,y),self.paths[-1]))
                    curr = (x,0)
                else:
                    (x,y) = (curr[0]-curr[1]/self.slope,0)
                    self.paths.append(Path(curr,(x,y),self.paths[-1]))
                    curr = (x,2)
            else:
                (x,y) = (0,curr[1]-self.slope*curr[0])
                if y <= 2 and y>=0:
                    self.paths.append(Path(curr,(x,y),self.paths[-1]))
                    curr = (x,2-y)
                    direction = not direction
                elif self.slope > 0:
                    (x,y) = (curr[0]-curr[1]/self.slope,0)
                    self.paths.append(Path(curr,(x,y),self.paths[-1]))
                    curr = (x,2)
                else:
                    (x,y) = (curr[0]+(2-curr[1])/self.slope,2)
                    self.paths.append(Path(curr,(x,y),self.paths[-1]))
                    curr = (x,0)
            if type(self.paths[-2])==type(self.paths[-1]):
                self.paths[-2].next = self.paths[-1]
        self.paths = self.paths[1:]

    def getc0(self):
        l = []
        i = 0
        #print(abs(i/self.slope))
        #print(abs(self.slope.denominator()))
        while abs(i/self.slope) <= abs(self.slope.denominator()):#Used to have denominator(), which was working before... Depends on if slope is an int or sage's data type
            if floor(i/self.slope)%2==0:
                l.append(i/self.slope-floor(i/self.slope))
            else:
                l.append(1+floor(i/self.slope)-i/self.slope)
            i+=2
        return(l)

    def getcpi(self):
        l = []
        i=1
        while abs(i/self.slope) <= abs(self.slope.denominator()):
            if floor(i/self.slope)%2==0:
                l.append(i/self.slope-floor(i/self.slope))
            else:
                l.append(1+floor(i/self.slope)-i/self.slope)
            #print(i)
            #print(l[-1])
            i+=2
        return(l)

class Path:
    def __init__(self,start,end,prev=None,nnext=None,flag=0,link=[],xhem=0):
        self.start = start
        self.end = end
        self.prev = prev
        self.next = nnext
        self.flag = flag #Flagged paths are colored differently in pictures
        self.slope = (self.end[1]-self.start[1])/(self.end[0]-self.start[0])

        #Only relevant for paths which represent nbd reps
        self.link = link #Lengths in the associated linkage for a nbd rep
        self.xhem = xhem #Which hemisphere x is sent to when b is sent to e^{\gamma k}i for 0<gamma<pi

    def reverse(self):
        (self.start,self.end) = (self.end,self.start)
        (self.prev,self.next) = (self.next,self.prev)

    def reslope(self):
        self.slope = (self.end[1]-self.start[1])/(self.end[0]-self.start[0])

    def perturb(self, d = 1/50):
        if self.prev.end == self.start:
            if self.prev.end[0]<self.prev.start[0]:
                pnudge = (self.start[0]+d,self.start[1]+d*self.prev.slope)
            else:
                pnudge = (self.start[0]-d,self.start[1]-d*self.prev.slope)
            self.start = pnudge
            self.prev.end = pnudge
            self.prev.reslope()
        if self.next.start == self.end:
            if self.next.end[0]>self.next.start[0]:
                nnudge = (self.end[0]+d,self.end[1]+d*self.next.slope)
            else:
                nnudge = (self.end[0]-d,self.end[1]-d*self.next.slope)
            self.end = nnudge
            self.next.start = nnudge
            self.next.reslope()

        self.reslope()


    def __str__(self):
        return str([self.start,self.end])

    def __repr__(self):
        return str([self.start,self.end])

class CRat(Tangle):
    #Cyclic tangle with rational subtangles
    def __init__(self,tangles):
        self.tangles = tangles
        self.paths = []
        self.gen()

    def gen(self):
        #Make BD reps
        slope = -1/sum([1/t.slope for t in self.tangles])
        T = RatTangle(slope)
        T.gen()
        T.breakall(1)
        self.paths = T.paths.copy()
        T.getcpi()
        #T.show()

        #NBD reps at c=0
        A = self.tangles[0].getc0()
        B = self.tangles[1].getc0()
        flag = 1
        for a in A:
            for b in B:
                #if a not in [0,1] and b not in [0,1]:
                if a>0 and a<1 and b>0 and b<1:
                    l = self.normx(a-b)
                    r = self.normx(a+b)
                    outs = self.getend((l,0))
                    oute = self.getstart((r,0))
                    ine = self.getstart((l,0))
                    ins = self.getend((r,0))

                    #Make sure that perturbing will always make them cross
                    if sign(outs.start[1]-outs.end[1]) == sign(oute.start[1]-oute.end[1]):
                        #Both ends are oriented the same (type Obar)
                        #xhem = 2*int((oute.end[1] > 1) ^ (a<b))-1 #Keeps track of the difference between the two arcs. Keeps track of which hemisphere x is sent to (in std representative)
                        xhem = 2*int(oute.end[1] > 1)-1 #Keeps track of the difference between the two arcs. Keeps track of which hemisphere x is sent to (in std representative)
                        pout = Path((l,0),(r,0),outs,oute,flag,[a,b],xhem)
                        outs.next = pout
                        oute.prev = pout

                        pin = Path((r,0),(l,0),ins,ine,flag+1,[a,b],-xhem)
                        ins.next = pin
                        ine.prev = pin
                        #print([a<b, 'O'])
                    else:
                        #Ends have opposite orientations (type S)
                        #xhem = 2*int((ins.end[1] > 1) ^ (a<b))-1 #Keeps track of the difference between the two arcs. Keeps track of which hemisphere x is sent to (in std representative)
                        xhem = 2*int(ins.end[1] > 1)-1 #Keeps track of the difference between the two arcs. Keeps track of which hemisphere x is sent to (in std representative)
                        pout = Path((l,0),(r,0),outs,ins,flag,[a,b],xhem)
                        outs.next = pout
                        ins.next = pout

                        pin = Path((r,0),(l,0),oute,ine,flag+1,[a,b],-xhem)
                        oute.prev = pin
                        ine.prev = pin
                        #print([a<b,'S'])
                    flag+=2

                    self.paths.append(pout)
                    self.paths.append(pin)

        #NBD reps at c=pi
        A = self.tangles[0].getcpi()
        B = self.tangles[1].getcpi()
        for a in A:
            for b in B:
                #if a not in [0,1] and b not in [0,1]:
                if a>0 and a<1 and b>0 and b<1:
                    l = self.normx(a-b)
                    r = self.normx(a+b)
                    #print([l,r])
                    outs = self.getend((l,1))
                    oute = self.getstart((r,1))
                    ine = self.getstart((l,1))
                    ins = self.getend((r,1))

                    if sign(outs.start[1]-outs.end[1]) == sign(oute.start[1]-oute.end[1]):
                        #Make sure that perturbing will always make them cross
                        #xhem = 2*int((oute.end[1] > 1) ^ (a<b))-1 #Keeps track of the difference between the two arcs. Keeps track of which hemisphere x is sent to (in std representative)
                        xhem = 2*int(oute.end[1] > 1)-1 #Keeps track of the difference between the two arcs. Keeps track of which hemisphere x is sent to (in std representative)
                        pout = Path((l,1),(r,1),outs,oute,flag,[a,b],xhem)
                        outs.next = pout
                        oute.prev = pout

                        pin = Path((r,1),(l,1),ins,ine,flag+1,[a,b],-xhem)
                        ins.next = pin
                        ine.prev = pin
                    else:
                        #xhem = 2*int((ins.end[1] > 1) ^ (a<b))-1 #Keeps track of the difference between the two arcs. Keeps track of which hemisphere x is sent to (in std representative)
                        xhem = 2*int(ins.end[1] > 1)-1 #Keeps track of the difference between the two arcs. Keeps track of which hemisphere x is sent to (in std representative)
                        pout = Path((l,1),(r,1),outs,ins,flag,[a,b],xhem)
                        outs.next = pout
                        ins.next = pout

                        pin = Path((r,1),(l,1),oute,ine,flag+1,[a,b],-xhem)
                        oute.prev = pin
                        ine.prev = pin
                    flag+=2

                    self.paths.append(pout)
                    self.paths.append(pin)

        self.split()
        #Reset orientations

    def perturb(self,d = 1/50):
        #Move NBD arcs vertically a little bit
        #NOTE Hacky way to get around the fact that I am changing the list as I iterate through it.
        i = 0
        while i < len(self.paths):
            p = self.paths[i]
            if p.flag:
                if p.start[1]==0 and p.end[1]==0:
                    (p1,p2) = self.subdivide(p,((p.start[0]+p.end[0])/2,0))
                    if p1.prev.end[1] == 2:
                        p1.start = (p1.start[0],2)
                        p1.end = (p1.end[0],2)
                    if p2.next.start[1] == 2:
                        p2.start = (p2.start[0],2)
                        p2.end = (p2.end[0],2)
                    p1.perturb(d)
                    p2.perturb(d)
                elif p.start[1]==p.end[1]:
                    p.perturb(d)
                    i+=1
                else:
                    i+=1
            else:
                i+=1

class CArb:
    #Cyclic tangle with arbitrary subtangles (or, unintentionally, arborescent)
    def __init__(self,tangles):
        self.tangles = tangles

    def gen(self):
        #Generates the character variety of the composite tangle based on the character varieties of the input tangles

        #Make sure that the parametrizationsa are the right ones
        t1 = self.tangles[0]
        t2 = self.tangles[1]

        sources = {}

        #TODO Could almost definitely optimize by modifying line-sweep
        #Generate the 
        for p1 in t1:
            [l1,h1] = sorted([p1.start,p1.end],key = lambda x: x[1])
            for p2 in t2:
                [l2,h2] = sorted([p2.start,p2.end],key = lambda x: x[1])
                l = max(l1[1],l2[1])
                h = min(h1[1],h2[1])
                if h <= l:
                    pass
                lx = (l-l1[1])/p1.slope+l1[0] + (l-l2[1])/p2.slope+l2[0]
                hx = h1[0] -(h1[1]-h)/p1.slope + h2[0] - (h2[1]-h)/p2.slope

                p = Path()

class IntersectionPoint:
    def __init__(self, coord, path1, path2, linkage = []):
        self.coord = coord
        self.p1 = path1
        self.p2 = path2
        self.grade = None
        if linkage == []:
            self.linkage = []
        else:
            self.linkage = linkage[linkage.index(min(linkage)):]+linkage[:linkage.index(min(linkage))] #Encodes the actual representation on C_3, useful to identify reps across different decompositions. In standard forn


    def __repr__(self):
        return str(self.coord)

class Bigon:
    def __init__(self, int1, int2, verts, intind = None):
        self.int1 = int1 #Source
        self.int2 = int2 #Target
        self.verts = verts
        self.twist = None
        self.lift = None
        self.word = None
        self.comps = None
        self.intind = intind
        self.endinlift = None #Coord of target int in the lift
        #self.comps = self.split()
        if len(self.verts) > 0:
            self.process()
            #self.orient()

    #TODO: Outdated?
    def split(self,verts = None):
        if verts == None:
            verts = self.verts
        vedge = None
        for n in range(len(verts)):
            (v1,v2) = (verts[n],verts[(n+1)%len(verts)])
            if equiv(v1,v2):
                if vedge == None:
                    vedge = v2
                elif v1[0]==vedge[0] or v1[1]==vedge[1]:
                    return [verts[verts.index(vedge):n+1]]+self.split(verts[:verts.index(vedge)]+verts[n+1:])
                else:
                    vedge = v2
        return [undup(verts)]

    def process(self):
        #Gets pi_1 of potential bigon on pillowcase. If pi_1 is trivial, gives a lift of bigon into R^2 and splits it into components on pillowcase to allow drawing.
        full = copy(self.verts)+[self.verts[0]]
        lift = [full[0]]
        adj = [0,0]
        flip = 1
        word = []
        splitind = [0]
        comps = [[]]
        for k in range(len(full)-1):
            comps[splitind[-1]].append(full[k])
            lift.append((flip*full[k][0]-adj[0],flip*full[k][1]-adj[1]))
            if k == self.intind:
            #if equiv(self.int2.coord, full[k]) and self.endinlift == None:
                #if self.endinlift != None:
                #    print([self.int2.coord,full[k],lift[-1]])
                self.endinlift = lift[-1]
            if full[k][0]==0 and full[k+1][0]==0 and full[k][1]+full[k+1][1]==2:
                #Left edge
                if word == [] or word[-1] != 0:
                    word.append(0)
                    comps.append([])
                    splitind.append(len(comps)-1)
                else:
                    word = word[:-1]
                    splitind = splitind[:-1]
                flip *= -1
                adj = [flip*full[k+1][0] - lift[-1][0],flip*full[k+1][1]-lift[-1][1]]
            elif full[k][0]==1 and full[k+1][0]==1 and full[k][1]+full[k+1][1]==2:
                #Right edge
                if word == [] or word[-1] != 1:
                    word.append(1)
                    comps.append([])
                    splitind.append(len(comps)-1)
                else:
                    word = word[:-1]
                    splitind = splitind[:-1]
                flip *= -1
                adj = [flip*full[k+1][0] - lift[-1][0],flip*full[k+1][1]-lift[-1][1]]
            elif full[k][0] == full[k+1][0] and set([full[k][1],full[k+1][1]])==set([0,2]):
                #Bottom edge
                if word == [] or word[-1] != 2:
                    word.append(2)
                    comps.append([])
                    splitind.append(len(comps)-1)
                else:
                    word = word[:-1]
                    splitind = splitind[:-1]
                adj = [flip*full[k+1][0] - lift[-1][0],flip*full[k+1][1]-lift[-1][1]]
            #if ((full[k][1]<1 and full[k+1][1]>=1) or (full[k][1]>=1 and full[k+1][1]<1)) and full[k][0]!=full[k+1][1]:
            elif ((full[k][1]<1 and full[k+1][1]>=1) or (full[k][1]>=1 and full[k+1][1]<1)):
                #Top (middle) edge
                #NOTE: This assumes no vertical lines which go all across the pillowcase. Could get around this by breaking at height 1 first and turning >= to ==
                if word == [] or word[-1] != 3:
                    word.append(3)
                else:
                    word = word[:-1]
            #else:
                #Add to lift
                #lift.append((flip*full[k][0]-adj[0],flip*full[k][1]-adj[1]))
        lift.append((flip*full[-1][0]-adj[0],flip*full[-1][1]-adj[1]))
        self.word = word
        self.lift = undup(lift)
        if word == []:
            self.comps = list(filter(lambda x: len(x)>2,map(undup,comps)))


    def orient(self):
        #Fix orientation by calculating winding number. In addition make sure that the angles at each intersection point are the correct sign
        """for c in self.comps:
            if self.int1.coord in c:
                ci = c
            if self.int2.coord in c:
                cj = c
        '''w = 0
        for n in range(len(ci)-2):
            w += findang(ci[n],ci[n+1],ci[n+2])
        w += findang(ci[-2],ci[-1],ci[0])
        ai = findang(ci[-1],ci[0],ci[1]) #Pretty sure int1 should always be the first vertex in its component NOTE: Should make sure this process preserves that property
        w += ai'''
        w = winding(ci)
        ww = winding(ci)

        print([N(w),N(ww)])

        if abs(w + 2*pi)<0.00001 or abs(ww + 2*pi)<0.00001: #Should be pm2pi, but there could be rounding error. Swap the inequality to switch conventions
            (self.int1,self.int2)=(self.int2,self.int1)

        #TODO: Check sign of ai and aj (once found)"""

        rev = False
        for c in self.comps:
            if len(c)<3:
                continue
            w = winding(c)
            if abs(w + 2*pi)<0.00001:
                #(self.int1,self.int2)=(self.int2,self.int1)
                #return
                rev = True
            elif abs(w - 2*pi)<0.00001:
                #return
                pass
            else:
                self.twist = True
                #print('Bad bigon')
                return
        if rev:
            (self.int1,self.int2)=(self.int2,self.int1)

    def z(self):
        #Calculates z as given in [HHK II] for computing gradings
        s=0
        for i in range(len(self.word))[::2]:
            l = self.word[i:i+2]
            if l in [[0,2],[2,1],[1,3],[3,0]]:
                s+=1
            elif l in [[2,0],[1,2],[3,1],[0,3]]:
                s-=1
            else:
                s+=2
        return s%4

    def mu(self,lfs=1+1/pi):
        #print([self.lift.index(self.endinlift),len(self.lift)])
        #print(self.lift)
        s = 0
        l = len(self.lift)
        if len(self.word)>0:
            l-=1
        for i in range(1,l):
            if self.lift[i] != self.endinlift:
                s1 = (self.lift[i][1]-self.lift[i-1][1])/(self.lift[i][0]-self.lift[i-1][0])
                s2 = (self.lift[(i+1)%len(self.lift)][1]-self.lift[i][1])/(self.lift[(i+1)%len(self.lift)][0]-self.lift[i][0])
                a = finddir(self.lift[i-1],self.lift[i],self.lift[(i+1)%len(self.lift)])

                if a < 0:
                    s -= t(s1,s2,lfs)
                elif a > 0:
                    s += t(s2,s1,lfs)
        return s%4

    def boundsimmerseddisk(self ,slope = N(1+1/pi)):
        #Implementation of Samuel Blanks algorithm (1976)
        #NOTE: Much faster if using approximation of irrational or rational number for slope
        #TODO: Could be substantially improved by using a sweep-line algorithm
        liftp = []
        for i in range(len(self.lift)):
            liftp.append(Path(self.lift[i],self.lift[(i+1)%len(self.lift)]))
        T = Tangle(liftp)
        #1. Split lift into regions
        # 1a. Find all self intersections and split paths there
        ints = []
        for n,p1 in enumerate(T.paths):
            for p2 in T.paths[n+2:]:
                i = pathint(p1,p2,atstart=False)
                if i:
                    ints.append(i)
        for i in ints:
            nps = T.subdivide(i.p1,i.coord)
            if nps != None:
                for j in ints:
                    if j.p1 == i.p1 and j != i:
                        if j.coord[0] >= min(nps[1].start[0],nps[1].end[0]) and j.coord[0] <= max(nps[1].start[0],nps[1].end[0]):
                            j.p1 = nps[1]
                        else:
                            j.p1 = nps[0]
                i.p1 = nps[1]
            nps = T.subdivide(i.p2,i.coord)
            if nps != None:
                for j in ints:
                    if j.p2 == i.p2 and j != i:
                        if j.coord[0] >= min(nps[1].start[0],nps[1].end[0]) and j.coord[0] <= max(nps[1].start[0],nps[1].end[0]):
                            j.p2 = nps[1]
                        else:
                            j.p2 = nps[0]
                i.p2 = nps[1]
        # 1b. 
        regions = []
        rev = deepcopy(T.paths)#Reversed orientation
        for rp in rev:
            rp.reverse()
        regp = T.paths + rev
        while len(regp)>0:
            reg = [regp[0]]
            while True:
                cands = list(filter(lambda x: x.start == reg[-1].end and x.end != reg[-1].start,regp))
                if len(cands)==1:
                    nextp = cands[0]
                else:
                    cands.sort(key = lambda x: findang(reg[-1].start, reg[-1].end, x.end))
                    nextp = cands[0]
                if nextp == reg[0]:
                    regions.append(reg)
                    for r in reg:
                        regp.remove(r)
                    break
                reg.append(nextp)

        if len(regions) == 2:
            return True
        if len(regions) == 3:
            return False
        # 1c. Remove the outside region. The winding number of the region around an interior point should all be the same EXCEPT for the outside region. Could also detect in the next loop and skip.

        #2. Find a point in each region and a line which goes from this point to outside the polygon
        lines = []
        for r in regions:
            #show(polygon([[yy.start[0],yy.start[1]] for yy in r]))
            #print([[N(yy.start[0]),N(yy.start[1])] for yy in r])
            bp = [(r[0].start[0]+r[0].end[0])/2,(r[0].start[1]+r[0].end[1])/2] #Point on the boundary of the region
            leftx = min([p[0] for p in self.lift])
            rightx = max([p[0] for p in self.lift])
            t1 = (leftx, bp[1]-(bp[0]-leftx)*slope)
            t2 = (rightx, bp[1]+(rightx-bp[0])*slope)
            trigger = Path(t1,t2)
            rints = []
            for p in r:
                ri = pathint(p,trigger)
                if ri:
                    rints.append(ri)
            rints.sort(key = lambda x: -x.coord[0])
            if finddir(rints[0].p1.start,rints[0].coord,t2) == 1:
                #Filters out the exterior region
                lines.append(Path(((rints[0].coord[0]+rints[1].coord[0])/2,(rints[0].coord[1]+rints[1].coord[1])/2),t2))
        #3. Find out where the polygon intersects the lines and turn this into a word
        wordgen = []
        for n,l in enumerate(lines):
            for m,p in enumerate(T.paths):
                i = pathint(l,p)
                if i:
                    nn = n+1 #Dont want 0 since 0 = -0
                    d = m+(i.coord[0]-p.start[0])/(p.end[0]-p.start[0])
                    if finddir(p.start,i.coord,l.end)<0:
                        #Seems like the wrong inequality, but the correct one for current orientation conventions
                        nn *= -1
                    wordgen.append([nn,d])
        word = [x[0] for x in sorted(wordgen,key = lambda x: x[1])]
        #4. Group the word
        return group(word)>0


    def getimg(self):
        im = Graphics()
        for c in self.comps:
            im += polygon(c,rgbcolor = 'red',alpha=0.5)
        return im


    def __repr__(self):
        return str(self.int1)+"->"+str(self.int2)
        #Currently direction is probably wrong

class Intersection:
    def __init__(self, t1, t2, distpath = [], iscomp = False):
        self.t1 = t1
        self.t2 = t2
        self.ips = []
        self.grips = [[],[],[],[]]
        self.bigons = []
        self.dist = None #Distinguished intersection
        self.iscomp = iscomp
        if not iscomp:
            for p1 in t1.paths:
                for p2 in t2.paths:
                    r = pathint(p1,p2)
                    if r == True:
                        print("Non-transverse!")
                    elif r != False:
                        self.ips.append(r)
                        if r.p1 in distpath or r.p2 in distpath:
                            self.dist = r
            for i in self.ips:
                #if not equiv(i.coord,i.p1.start):
                nps = self.t1.subdivide(i.p1,i.coord)
                if nps != None:
                    for j in self.ips:
                        if j.p1 == i.p1 and j != i:
                            if j.coord[0] >= min(nps[1].start[0],nps[1].end[0]) and j.coord[0] <= max(nps[1].start[0],nps[1].end[0]):
                                j.p1 = nps[1]
                            else:
                                j.p1 = nps[0]
                    i.p1 = nps[1]
                nps = self.t2.subdivide(i.p2,i.coord)
                if nps != None:
                    for j in self.ips:
                        if j.p2 == i.p2 and j != i:
                            if j.coord[0] >= min(nps[1].start[0],nps[1].end[0]) and j.coord[0] <= max(nps[1].start[0],nps[1].end[0]):
                                j.p2 = nps[1]
                            else:
                                j.p2 = nps[0]

                    i.p2 = nps[1]
            #CHECK
            for i in self.ips:
                if not equiv(i.coord,i.p1.start) or not equiv(i.coord,i.p2.start):
                    print('BAD split')
                    print([i,i.p1,i.p2])


            self.comps = self.split()

        else:
            self.comps = [self]

    def show(self):
        im = self.t1.getimg()+self.t2.getimg(colup = 'red',coldown='red')
        for i in self.ips:
            im += circle(i.coord,0.0125,rgbcolor='black',fill=True)
        show(im,aspect_ratio=1,figsize=12,tick_formatter=[formatter,formatter], ticks=[1/2,1/2], axes_labels=[r'$\gamma$',r'$\theta$'])

    def showcomps(self):
        ga = []
        for c in self.comps:
            im = c.t1.getimg()+c.t2.getimg(colup = 'red',coldown='red')
            for i in c.ips:
                im += circle(i.coord,0.0125,rgbcolor='black',fill=True)
            ga.append(im)
        show(graphics_array(ga),aspect_ratio=1,figsize=12,tick_formatter=[formatter,formatter], ticks=[1/2,1/2], axes_labels=[r'$\gamma$',r'$\theta$'])

    def showgr(self):
        im = self.t1.getimg()+self.t2.getimg(colup = 'red',coldown='red')
        for i in self.grips[3]:
            im += circle(i.coord,0.025,rgbcolor=(.4,.4,1),fill=True)
        for i in self.grips[2]:
            im += circle(i.coord,0.019,rgbcolor=(0,1,0),fill=True)
        for i in self.grips[1]:
            im += circle(i.coord,0.013,rgbcolor=(1,0,0),fill=True)
        for i in self.grips[0]:
            im += circle(i.coord,0.007,rgbcolor='black',fill=True)
        show(im,aspect_ratio=1,figsize=12,tick_formatter=[formatter,formatter], ticks=[1/2,1/2], axes_labels=[r'$\gamma$',r'$\theta$'])

    def showbigon(self,b):
        im = self.t1.getimg()+self.t2.getimg(colup = 'red',coldown='red') + b.getimg()
        for i in [b.int1,b.int2]:
            if i.grade == 3:
                im += circle(i.coord,0.025,rgbcolor=(.4,.4,1),fill=True)
            if i.grade == 2:
                im += circle(i.coord,0.019,rgbcolor=(0,1,0),fill=True)
            if i.grade == 1:
                im += circle(i.coord,0.013,rgbcolor=(1,0,0),fill=True)
            if i.grade == 0:
                im += circle(i.coord,0.007,rgbcolor='black',fill=True)
        lim = Graphics() #Lift image
        #lim.axes_range(floor(min([x[0] for x in b.lift])), ceil(max([x[0] for x in b.lift])), floor(min([x[1] for x in b.lift])), ceil(max([x[0] for x in b.lift])))
        lim += point((floor(min([x[0] for x in b.lift])), floor(min([x[1] for x in b.lift]))), color='white')
        lim += point((ceil(max([x[0] for x in b.lift])), ceil(max([x[1] for x in b.lift]))), color='white')
        lim += polygon([[x[0],x[1]] for x in b.lift])
        show(graphics_array([im,lim]),aspect_ratio=1,figsize=12,tick_formatter=[formatter,formatter], ticks=[1/2,1/2], axes_labels=[r'$\gamma$',r'$\theta$'])
        #show(im,aspect_ratio=1,figsize=12,tick_formatter=[formatter,formatter], ticks=[1/2,1/2], axes_labels=[r'$\gamma$',r'$\theta$'])


    def showchains(self):
        p = point((0,0),rgbcolor=(1,1,1),axes=False)
        for n,i in enumerate(self.grips[0]):
            n+=1
            p+=point((0,n),rgbcolor=(1-Integer(len(i.linkage)==1),0,Integer(len(i.linkage)==1)))
            for b in self.bigons:
                if b.int1 == i:
                    p+=line([(0,n),(self.grips[1].index(b.int2)+1,0)])
        for n,i in enumerate(self.grips[1]):
            n+=1
            p+=point((n,0),rgbcolor=(1-Integer(len(i.linkage)==1),0,Integer(len(i.linkage)==1)))
            for b in self.bigons:
                if b.int1 == i:
                    p+=line([(n,0),(0,-self.grips[2].index(b.int2)-1)])
        for n,i in enumerate(self.grips[2]):
            n+=1
            p+=point((0,-n),rgbcolor=(1-Integer(len(i.linkage)==1),0,Integer(len(i.linkage)==1)))
            for b in self.bigons:
                if b.int1 == i:
                    p+=line([(0,-n),(-self.grips[3].index(b.int2)-1,0)])
        for n,i in enumerate(self.grips[3]):
            n+=1
            p+=point((-n,0),rgbcolor=(1-Integer(len(i.linkage)==1),0,Integer(len(i.linkage)==1)))
            for b in self.bigons:
                if b.int1 == i:
                    p+=line([(-n,0),(0,self.grips[0].index(b.int2)+1)])
        show(p)

    def split(self):
        #Splits Tangle 1 into multiple intersections
        comps = self.t1.split()
        if len(comps)==1:
            self.iscomp = True
            return [self]
        else:
            ints = []
            for n,c in enumerate(comps):
                i = Intersection(c,self.t2,iscomp = True)
                for j in self.ips:
                    if j.p1 in c.paths:
                        i.ips.append(j)
                #for b in self.bigons:
                #    if b.int1 in i.ips:
                #        i.bigons.append(b)
                if n==0:
                    i.dist = self.dist
                ints.append(i)
            return ints

    def findbigons(self):
        for p in self.ips:
            for q in self.grips[(p.grade+1)%4]:
        #for i,p in enumerate(self.ips):
            #for q in self.ips[i+1:]:
                p1f = findpath(p.p1,q.p1)
                p1b = findpath(p.p1,q.p1,True)
                v1 = []
                if p1f != None:
                    v1.append([y for x in p1f for y in [x.start,x.end]]+[q.coord])
                if p1b != None:
                    v1.append([p.coord]+[y for x in p1b for y in [x.end,x.start]])

                p2f = findpath(p.p2,q.p2)
                p2b = findpath(p.p2,q.p2,True)
                v2  = []
                if p2f != None:
                    v2.append([y for x in p2f for y in [x.start,x.end]]+[q.coord])
                if p2b != None:
                    v2.append([p.coord]+[y for x in p2b for y in [x.end,x.start]])
                for a in v1:
                    for b in v2:
                        b = Bigon(p,q,undup(a+b[::-1]))
                        if len(b.word)==0:
                            if b.boundsimmerseddisk():
                                self.bigons.append(b)

    def homology(self):
        gens = copy(self.ips)
        arrows = copy(self.bigons)
        while len(arrows)>0:
            #Delete arrow and add zig-zags
            a = arrows[0]
            zz1 = filter(lambda x: x.int2 == a.int2 and x != a,arrows)
            zz2 = filter(lambda x: x.int1 == a.int1 and x != a,arrows)
            #print(list(zz1))
            #print(list(zz2))
            for za in zz1:
                for zb in zz2:
                    if za == zb:
                        print('Bad chain complex')
                        return
                    arrows.append(Bigon(za.int1,zb.int2,[]))
            arrows = list(filter(lambda x: x.int1 not in [a.int1,a.int2] and x.int2 not in [a.int1,a.int2],arrows))
            if a.int1==a.int2: #For debugging
                print("Bad chain complex")
                return []
            gens.remove(a.int1)
            gens.remove(a.int2)
        grhom = [len(list(filter(lambda x: x.grade == i,gens))) for i in range(4)]
        return grhom

    def grgens(self):
        return [len(x) for x in self.grips]

    def grade(self,lfslope=1+1/pi):
        if self.iscomp:
            self.grcorr = "?"
            if self.dist != None:
                d = self.dist
                self.grcorr = "sigma"
            else:
                d = self.ips[0]
            for p in self.t1.paths+self.t2.paths:
                if p.slope == lfslope:
                    print("Bad linefield slope")
                    return
            
            debug = []
            for i in self.ips:
                if d == i:
                    i.grade = 0 #+correction term
                    gr = 0
                else:
                    #Calculate relative grading gr(d,i)
                    #Find path in Lagrangian1
                    p1 = findpath(d.p1,i.p1)
                    #print([d,i,p1])
                    if p1 != None:
                        v1 = [y for x in p1 for y in [x.start,x.end]]+[i.coord]
                    else:
                        p1 = findpath(d.p1,i.p1,True)
                        v1 = [d.coord]+[y for x in p1 for y in [x.end,x.start]]

                    s1d = p1[0].slope
                    s1i = p1[-1].slope

                    p2 = findpath(d.p2,i.p2)
                    if p2 != None:
                        v2 = [y for x in p2 for y in [x.start,x.end]]+[i.coord]
                    else:
                        p2 = findpath(d.p2,i.p2,True)
                        #v2 = [d.coord]+sum([[x.end,x.start] for x in p2])
                        v2 = [d.coord]+[y for x in p2 for y in [x.end,x.start]]
                    s2d = p2[0].slope
                    s2i = p2[-1].slope

                    ts = t(s1d,s2d,lfslope)-t(s1i,s2i,lfslope)
                    #Make bigon, get z, use lift to get 
                    fullv = undup(v1)
                    endind = len(fullv)-1
                    if v1[0]==v1[-1]:
                        #Since undup cuts off the last part of the path if it goes all the way around
                        endind += 1
                    fullv = undup(fullv+v2[::-1])
                    #print([endind,fullv])
                    b = Bigon(d,i,fullv,intind = endind)
                    #b.process()
                    z = b.z()
                    mu = b.mu(lfslope)
                    gr = (ts+mu+z)%4 #TODO:Signs?
                    #print([ts,mu,z,b.word])
                    #print(p1+p2)
                    #print(p1)
                    #print(v1)
                    i.grade = gr
                    debug.append(line(b.lift,axes=False)+point(b.lift[0],rgbcolor=(1,0,0))+point(b.endinlift,rgbcolor=(0,1,0)))
                self.grips[gr].append(i)
            #show(graphics_array(debug[:5]))
        else:
            for c in self.comps:
                c.grade(lfslope)
            for i in self.ips:
                self.grips[i.grade].append(i)

    def identify(self,ref):
        #Promote the relative grading given a different splitting which already has an absolute grading. Also reorder grips so they can be compared.
        newgrips = [[None]*len(ref.grips[0]),[None]*len(ref.grips[1]),[None]*len(ref.grips[2]),[None]*len(ref.grips[3])]
        for i in self.comps[0].ips:
            if i == self.dist:
                match = ref.dist
            else:
                #print([i,i.linkage])
                #print(list(filter(lambda x: x.linkage == i.linkage, ref.ips)))
                matches = list(filter(lambda x: x.linkage == i.linkage and x.grade == i.grade, ref.ips))
                if len(matches) != 1:
                    print('Non-unique match')
                match = matches[0]
            newgrips[i.grade][ref.grips[i.grade].index(match)] = i
        for c in self.comps[1:]:
            gradeshift = None
            for i in c.ips:
                print([i,i.grade])
                matches = list(filter(lambda x: x.linkage == i.linkage, ref.ips))
                print([[x,x.grade] for x in matches])
                if len(matches) == 1:
                    gradeshift = (matches[0].grade-i.grade)%4
                    #break
                elif len(matches) == 2:
                    mgrades = [x.grade for x in matches]
                    selfmatches = list(filter(lambda x: x.linkage == i.linkage, self.ips))
                    smgrades = [x.grade for x in selfmatches]
                    for n in range(4):
                        if set(mgrades) == set([(x+n)%4 for x in smgrades]):
                            gradeshift = n
                    if gradeshift == None:
                        print('Hmm')
                    #break
                else:
                    print('Hmmm')
                print(gradeshift)
            for i in c.ips:
                print(gradeshift)
                i.grade = (i.grade+gradeshift)%4
                matches = list(filter(lambda x: x.linkage == i.linkage and x.grade == i.grade, ref.ips))
                if len(matches) != 1:
                    print('Non-unique match')
                match = matches[0]
                newgrips[i.grade][ref.grips[i.grade].index(match)] = i
        self.grips = newgrips




def equiv(p1,p2):
    #Determines if two points on the pillocase are equivalent
    if p1==p2:
        return True
    if p1[0]==p2[0] and p1[0] in [0,1] and p1[1]+p2[1]==2:
        return True
    if p1[0]==p2[0] and set([p1[1],p2[1]])==set([0,2]):
        return True
    return False

def findang(p1,p2,p3):
    #Find the signed angle between three points
    v = [p2[0]-p1[0],p2[1]-p1[1]]
    w = [p3[0]-p2[0],p3[1]-p2[1]]
    ang = acos((v[0]*w[0]+v[1]*w[1])/sqrt((v[0]**2+v[1]**2)*(w[0]**2+w[1]**2)))
    d = finddir(p1,p2,p3)
    if d<0:
        ang *= -1
    return ang

def finddir(p1,p2,p3):
    #Special case of findang, but since that was the time bottleneck, this should speed things up
    k = (p2[0]-p1[0])*(p3[1]-p2[1])-(p2[1]-p1[1])*(p3[0]-p2[0])
    return sign(k)


def winding(v):
    w = 0
    for n in range(len(v)-2):
        w += findang(v[n],v[n+1],v[n+2])
    w += findang(v[-2],v[-1],v[0])
    ai = findang(v[-1],v[0],v[1])
    w += ai
    return w

def t(s1,s2,s):
    #Computes the triple index as defined in [HHK II]
    if s2<s1:
        if s2<s and s<s1:
            return 1
        else:
            return 0
    if s1<s2:
        if s2<s or s<s1:
            return 1
        else:
            return 0

def mu(p,lfs):
    pass

def findpath(source, target, back = False):
    p = source
    c = [p]
    while True:
        if back:
            p = p.prev
        else:
            p = p.next
        if p == None or (p == source and source != target):
            return None
        c.append(p)
        if p == target:
            if back:
                return c[1:]
            else:
                return c[:-1]

def pathint(p1,p2,atstart = True):
    #Check a loop for self-intersections. If atstart, normalize so all intersections on endpoints occur at the start
    if p1.slope != p2.slope:
    #Ignore parallel for now. Will probably need to fix later in case things are not transverse
        s = (p1.start[1]+p1.slope*(p2.start[0]-p1.start[0])-p2.start[1])/(p2.slope-p1.slope)
        i = (p2.start[0]+s,p2.start[1]+p2.slope*s)
        if i[0]>=min(p1.start[0],p1.end[0]) and i[0]<=max(p1.start[0],p1.end[0]) and i[0]>=min(p2.start[0],p2.end[0]) and i[0]<=max(p2.start[0],p2.end[0]):
            if atstart:
                if i == p1.end:
                    #Makes it so intersections are always at the beginning of paths
                    p1 = p1.next
                    i = p1.start
                if i == p2.end:
                    p2 = p2.next
                    i = p2.start
            if p1.xhem == 1:
                linkage = [p1.link[0], p1.link[1], i[0]]
            elif p1.xhem == -1:
                linkage = [p1.link[1], p1.link[0], i[0]]
            else:
                #For bd reps, c must match. Will be difficult to extract actual linkage
                if i[1]<1:
                    linkage = [i[1]] #For bd reps, c must match. Will be difficult to extract actual linkage
                else:
                    linkage = [2-i[1]]
            ip = IntersectionPoint(i,p1,p2,linkage)
            return ip
            '''elif i == p2.end and i != p2.next.start:
            #Edge case where two paths instersect at an edge of the pillowcase in a way where it wouldn't be normally counted
            print(["HERE",p1,p2,i])
            ip = IntersectionPoint(p2.next.start,p1,p2.next)
            return ip
            elif i == p1.end and i != p1.next.start:
            print(["HERE",p1,p2,i])
            ip = IntersectionPoint(p1.next.start,p1.next,p2)
            return ip
            #NOTE: Could double count if two paths cross transversely on an edge of the pillowcase with 'opposite' orientations'''
        else:
            return False
    else:
        #if (p2.start[1]-p1.start[1])/(p2.start[0]-p1.start[0])==p1.slope and (dist(p1.start,p2.start)+dist(p2.start,p1.end) <= dist(p1.start,p1.end) or dist(p1.start,p2.end)+dist(p2.end,p1.end)<=dist(p1.start,p1.end)):
        if dist(p1.start,p2.start)+dist(p2.start,p1.end) <= dist(p1.start,p1.end) or dist(p1.start,p2.end)+dist(p2.end,p1.end)<=dist(p1.start,p1.end):
            #Non-transverse intersection
            return True
        else:
            return False

def dist(v1,v2):
    return sqrt((v2[0]-v1[0])**2 + (v2[1]-v1[1])**2)

def undup(l):
    #Remove consecutive elements of a list which are equivalent
    t = []
    for e in l:
        if t == [] or e != t[-1]:
            t.append(e)
    if t[0] == t[-1]:
        t = t[:-1]
    return t

def group(l):
    #Implementation of grouping as described in [Blank]. Returns the number of groupings of a word. If the number is positive, then there is a way to extend the immersion of the circle which generated the word to an immersion of a disk whose boundary is the circle.
    #l is a list of numbers where negatives represent the inverses
    if len(l) == 0 or min(l)>0:
        #l is a positive word.
        return 1
    s=0
    for i in range(len(l)):
        for j in range(i+1,len(l)):
            if l[i] == -l[j]:
                if len(l[i+1:j]) == 0 or min(l[i+1:j])>0:
                    s += group(l[:i]+l[j+1:])
    return s


def chaincx(p,q,r, ear = 1/200):
    #Genreate chain complex. Assume if any are odd, it is the first entry
    T = RatTangle(q)
    S = RatTangle(r)
    C = CRat([T,S])
    E = RatTangle(p)
    dp = E.earring(ear)
    I = Intersection(C,E,distpath = dp)
    I.grade()
    return I.grgens()

def monthom(p,q,r, v = False, b = False, ear = 1/30, pert = 1/50):
    #Finds the Lagrangian Floer homology associated to a Montesinos knots
    T = RatTangle(p)
    S = RatTangle(q)
    C = CRat([T,S])
    if pert != 0:
        C.perturb(pert)
    E = RatTangle(r)
    dp=E.earring(ear)
    I = Intersection(C,E,distpath=dp)
    #return I


    print('Gens: '+str(len(I.ips)))
    I.grade()
    print('Gradings: '+str(I.grgens()))

    I.findbigons()
    
    bc = len(I.bigons)

    s = I.homology()

    print('Bigons: '+str(bc))
    print('Hom: '+str(s))
    print('Total Rank: '+str(sum(s)))
    print(' ')
    if v:
        I.showcomps()
        I.showgr()
        I.showchains()
        if b:
            for bg in I.bigons:
                I.showbigon(bg)
    return I

def rotate(p,q,r,pert=0,ear=1/1000,v=False,b=False,compare=False):
    #Computes the Lagrangian Floer homology of a Montesinos knot with 3 different decompositions and compares the resulting homologies
    I = monthom(p,q,r,pert=pert,ear=ear,v=v,b=b)
    J = monthom(q,r,p,pert=pert,ear=ear,v=v,b=b)
    K = monthom(r,p,q,pert=pert,ear=ear,v=v,b=b)
    print()
    if compare:
        print(sum(I.homology())==sum(J.homology()) and sum(J.homology())==sum(K.homology()))
        print()

def monttosnappylink(mont):
    #For some reason, sum doesn't work
    K = snappy.RationalTangle(0)
    for T in map(lambda x: snappy.RationalTangle(x.numerator(),x.denominator()), mont):
        K += T
    return K.numerator_closure()

def montsig(mont):
    k = monttosnappylink(mont)
    return k.signature()

def pretzsig(pqr):
    mont = map(lambda x: 1/x, pqr)
    return montsig(mont)

if __name__ == '__main__':
    pass
