import numpy as np
import math
def read_file(file_path):
    '''
    读取一个txt文件
    返回一个类，类的key就是txt中带有#的行的第一个内容
        比如#BUS,k，则key是'BUS'
    这一带#的行下面的所有内容都将被以list的形式存入content['BUS']中，
    '''
    content = {}
    key = ''

    single = False
    content_type = str
    with open(file_path, 'r') as f:
        lines = f.readlines()
    for line in lines:
        if line[0]=='#':
            line_split = line[1:-1].replace(' ','').split(',')
            key = line_split[0]
            if 'single' in line_split:
                single = True
                content[key] = None
            else:
                single = False
                content[key] = []
            if 'int' in line_split:
                content_type = int
            elif 'float' in line_split:
                content_type = float
            else:
                content_type = str
        elif line.replace('\n', '').replace(' ', ''):
            if single:
                if content_type == str:
                    content[key] = line.replace('\n', '').replace(' ', '')
                else:
                    content[key] = content_type(line.replace('\n', '').replace(' ', ''))
            else:
                line_content = line.replace('\n', '').replace(' ', '').split(',')
                if content_type == str:
                    content[key].append(line_content)
                else:
                    content[key].append([content_type(x) for x in line_content])
    return content

content = read_file("IEEE39node.txt")



class System:

    def __init__(self,content):
        self.load_num = 1
        self.version_str = content['VER']  # 版本号，字符串         1
        self.system_name = content['SYS']  # 系统名，字符串         IEEE39bussystem
        self.node_num, self.line_num, self.gen_num = content['INF'][0]  # 记录系统负荷、线路、发电机、变压器数量 39 46 10
        self.baseMVA = content['BAS']  # 基准功率值     100MVA
        self.line = content['LIN']
        self.bus = content['BUS']
        self.Gen = content['GEN']
        self.bus_type1 = 0        # self.bus_type123分别存放三类节点的个数
        self.bus_type2 = 0
        self.bus_type3 = 0
        self.node_type3 = 0
        for bus_info in self.bus:
            inf0,inf1,inf2,inf3,inf4,inf5 = bus_info
            node_num = int(inf0)
            type = int(inf1)
            if type ==1:
                self.bus_type1 +=1
            elif type ==2 :
                self.bus_type2 +=1
            elif type ==3 :
                self.bus_type3 +=1
                self.node_type3 = node_num      # self.node_type3 表示参考节点是第几个
            else:
                pass
        self.E0 = np.zeros([self.node_num, 1])
        self.F0 = np.zeros([self.node_num, 1])
        self.a = np.zeros([self.node_num,1])     # self.a,b用于计算功率的误差值
        self.b = np.zeros([self.node_num,1])
        self.P = np.zeros([self.bus_type1 + self.bus_type2, 1])
        self.P_delta = np.zeros([self.bus_type1 + self.bus_type2, 1])
        self.Q = np.zeros([self.bus_type1, 1])
        self.Q_delta = np.zeros([self.bus_type1, 1])
        self.U_2 = np.zeros([self.bus_type2, 1])
        self.U_2_delta = np.zeros([self.bus_type2, 1])


    # fp fq fu 用于进行节点坐标的变换，解决参考节点混在第一，第二类节点中，计算时标号的变化
    def fp(self,p):   #实现坐标转换   p的范围应该在 range(0,38) 即n-1个点
        if p < self.node_type3-1:
            return p
        else:
            return p+1

    def fq(self,q):  #实现q的坐标转化，q的范围应在 range(0,28)  即 n-r-1个点
        if self.node_type3 <= self.bus_type1:
            if q<self.node_type3-1:
                return q
            else:
                return q+1
        else:
            return q

    def fu(self,u):
        if self.node_type3<= self.bus_type1+1:
            return u+self.bus_type1 + 1
        else:
            if u < self.node_type3 - self.bus_type1 - 1:
                return u + self.bus_type1
            else:
                return u + self.bus_type1 + 1

    # YBus 用于获得系统的节点导纳矩阵，结果用self.G 和 self.B表示实部和虚部
    def YBus(self):
        self.G = np.zeros([self.node_num, self.node_num])
        self.B = np.zeros([self.node_num, self.node_num])
        for line_info in self.line:
            inf0, inf1, inf2, inf3, inf4, inf5 = line_info
            bus0, bus1, r, x, b, tap = eval(inf0), eval(inf1), eval(inf2), eval(inf3), eval(inf4), eval(inf5)
            if tap == 0:   # 线路无变压器
                self.G[bus0 - 1, bus1 - 1] = -r / (r * r + x * x)
                self.G[bus1 - 1, bus0 - 1] = -r / (r * r + x * x)
                self.B[bus0 - 1, bus1 - 1] = x / (r * r + x * x)
                self.B[bus1 - 1, bus0 - 1] = x / (r * r + x * x)  # 非对角元素
                self.G[bus0 - 1, bus0 - 1] += r / (r * r + x * x)
                self.G[bus1 - 1, bus1 - 1] += r / (r * r + x * x)
                self.B[bus0 - 1, bus0 - 1] += (b / 2 - x / (r * r + x * x))
                self.B[bus1 - 1, bus1 - 1] += (b / 2 - x / (r * r + x * x))
            else:         # 线路有变压器
                self.G[bus0 - 1, bus1 - 1] = -r / ((r * r + x * x) * tap)
                self.G[bus1 - 1, bus0 - 1] = -r / ((r * r + x * x) * tap)
                self.B[bus0 - 1, bus1 - 1] = x / ((r * r + x * x) * tap)
                self.B[bus1 - 1, bus0 - 1] = x / ((r * r + x * x) * tap)
                self.G[bus0 - 1, bus0 - 1] += r / ((r * r + x * x) * tap * tap)
                self.G[bus1 - 1, bus1 - 1] += r / (r * r + x * x)
                self.B[bus0 - 1, bus0 - 1] += (b / 2 - x / (r * r + x * x)) * (1 / (tap * tap))
                self.B[bus1 - 1, bus1 - 1] += b / 2 - x / (r * r + x * x)

    # 获得 设置各个节点的电压初值，非参考节点设置ei=1,fi=0
    def Get_EF0(self):                    #为e,f设定初始值

        # t = 0
        for bus_info in self.bus:
            inf0, inf1, inf2, inf3, inf4, inf5 = bus_info
            node, node_type, p, q, V_base, V_init = eval(inf0), eval(inf1), eval(inf2), eval(inf3), eval(inf4), eval(
                inf5)
            if node_type != 3:
                self.E0[node-1, 0] = 1
                self.F0[node-1, 0] = 0
                # t = t + 1
            else:
                self.E0[node-1, 0] = V_init
                self.F0[node-1, 0] = 0
                # t = t + 1

    # 获得用于计算功率误差的ai,bi
    def Get_AB(self):
        for i in range(self.node_num):
            for j in range(self.node_num):
                self.a[i,0] += (self.G[i,j]*self.E0[j,0] - self.B[i,j]*self.F0[j,0])
                self.b[i,0] += (self.G[i,j]*self.F0[j,0] + self.B[i,j]*self.E0[j,0])

    # 获得功率的输入值
    def Get_PQU0ic(self):
        pt = 0
        qt = 0
        ut = 0
        flag = 0
        for bus_info in self.bus:
            inf0, inf1, inf2, inf3, inf4, inf5 = bus_info
            node, node_type, p, q, V_base, V_init = eval(inf0), eval(inf1), eval(inf2), eval(inf3), eval(inf4), eval(
                inf5)
            if node_type==1 :
                self.P[pt,0] = -1*p/self.baseMVA
                self.Q[qt,0] = -1*q/self.baseMVA
                pt += 1
                qt += 1
            elif node_type==2 :
                ####
                # print(type(p))
                ####
                self.P[pt,0] = -1*p/self.baseMVA + float(float(self.Gen[ut+flag][1])/self.baseMVA)     # 还应考虑发电机
                pt += 1
                self.U_2[ut,0] = V_init
                # self.
                ut += 1
            else:
                flag = 1

    # 计算功率误差
    def DPQC(self):
        for i in range(self.bus_type1+self.bus_type2):
            xp = self.fp(i)
            self.P_delta[i,0] = self.P[i,0] - (self.E0[xp,0]*self.a[xp,0] + self.F0[xp,0]*self.b[xp,0])

        for i in range(self.bus_type1):
            xq = self.fq(i)
            self.Q_delta[i,0] = self.Q[i,0] - (self.F0[xq,0]*self.a[xq,0] - self.E0[xq,0]*self.b[xq,0])

        for i in range(self.bus_type2):
            xu = self.fu(i)
            # print(xu)
            # print(f'{self.U_2[i,0]}')
            self.U_2_delta[i,0] = self.U_2[i,0]*self.U_2[i,0] - (self.E0[xu,0]*self.E0[xu,0] + self.F0[xu,0]*self.F0[xu,0])


        #计算功率误差

    # 计算雅克比矩阵
    def JMCC(self):
        self.J = np.zeros([2*(self.bus_type1+self.bus_type2),2*(self.bus_type1+self.bus_type2)])
        #先生成H矩阵
        for i in range(self.bus_type1+self.bus_type2):
            for j in range(self.bus_type1+self.bus_type2):
                if self.fp(i)==self.fp(j):
                    x = self.fp(i)
                    # 生成分块矩阵H
                    self.J[i,i] = -1*self.a[x,0]-(self.G[x,x]*self.E0[x,0]+self.B[x,x]*self.F0[x,0])
                    # 生成分块矩阵N
                    self.J[i,i+self.bus_type1+self.bus_type2] = -1*self.b[x,0] + (self.B[x,x]*self.E0[x,0] - self.G[x,x]*self.F0[x,0])
                else:
                    xi = self.fp(i)
                    xj = self.fp(j)
                    # H矩阵的非对角元
                    self.J[i,j] = -1*(self.G[xi,xj]*self.E0[xi,0]+self.B[xi,xj]*self.F0[xi,0])
                    # N矩阵的非对角元
                    self.J[i,j+self.bus_type1+self.bus_type2] = self.B[xi,xj]*self.E0[xi,0]-self.G[xi,xj]*self.F0[xi,0]
        for i in range(self.bus_type1):        # n-1-r
            for j in range(self.bus_type1+self.bus_type2):   # n-1
                if self.fq(i) == self.fp(j):
                    #生成M矩阵的对角元素，非对角元素可由N矩阵得到
                    x = self.fq(i)
                    self.J[i+self.bus_type1+self.bus_type2,j] = self.b[x,0]+(self.B[x,x]*self.E0[x,0]-self.G[x,x]*self.F0[x,0])
                    # 生成L矩阵的对角元素，非对角元素可由H矩阵得到
                    self.J[i+self.bus_type1+self.bus_type2,j+self.bus_type2+self.bus_type1] = -1*self.a[x,0]+(self.G[x,x]*self.E0[x,0]+self.B[x,x]*self.F0[x,0])
                else:
                    # M矩阵的非对角元等于N矩阵同样位置的非对角元
                    self.J[i+self.bus_type1+self.bus_type2,j] = self.J[i,j+self.bus_type1+self.bus_type2]
                    # L矩阵的非对角元等于H矩阵同样位置的非对角元的相反数
                    self.J[i+self.bus_type1+self.bus_type2,j+self.bus_type2+self.bus_type1] = -1*self.J[i,j]

        for i in range(self.bus_type2):       # r
            for j in range(self.bus_type2+self.bus_type1):   # n-1
                if self.fu(i)==self.fp(j):
                    x = self.fu(i)
                    self.J[i+2*self.bus_type1+self.bus_type2,j] = -2*self.E0[x,0]
                    self.J[i+2*self.bus_type1+self.bus_type2,j+self.bus_type1+self.bus_type2] = -2*self.F0[x,0]
                else:
                    pass
        #形成雅克比矩阵

    # 求解方程，并且更新各个节点的电压值
    def SEVC(self):           # 解修正方程并且修改E0 F0矩阵
        self.E_delta = np.zeros([self.bus_type1+self.bus_type2,1])
        self.F_delta = np.zeros([self.bus_type1+self.bus_type2,1])
        J_1 = np.linalg.inv(self.J)
        t1 = np.append(self.P_delta,self.Q_delta,axis=0)
        t = np.append(t1,self.U_2_delta,axis=0)
        R = np.dot(J_1,t)
        self.E_delta,self.F_delta = np.vsplit(R,2)
        t = 0
        for i in range(self.node_num):
            if i != self.node_type3-1:
                self.E0[i,0] = self.E0[i,0] - self.E_delta[t,0]
                self.F0[i,0] = self.F0[i,0] - self.F_delta[t,0]
                t = t + 1
            else:
                pass



        # self.E0 = np.add(self.E0,-self.E_delta)
        # self.F0 = np.add(self.F0,-self.F_delta)
        #解修正方程

#创建对象，并获得系统的节点导纳矩阵，并为参考节点以外的节点电压赋初值

sys39 = System(content)
sys39.YBus()
sys39.Get_EF0()

# 以下五个函数可以完成一次叠代
sys39.Get_AB()
sys39.Get_PQU0ic()
sys39.DPQC()
sys39.JMCC()
sys39.SEVC()   # 算出并修改生成新的E0 F0矩阵

