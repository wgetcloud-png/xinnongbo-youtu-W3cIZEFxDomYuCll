
#### 提纲


1. 背景介绍
2. 三角方程组
3. Gauss消去法
4. 附录


## 一、背景介绍


### 1\.1 线性方程组的相关概念


线性方程组在解决现实实际问题中直接产生，最小二乘数据拟合、微分方程边值问题和初边值问题的数值解产生了大量的线性方程组。
线性方程组系数矩阵的类型分别有


1. **稠密型**(dense)：几乎所有元素都是非零的
2. **稀疏型**(sparse)：有大量零元素
3. **带状的**(banded)
4. **三角状**(triangular)
5. **块状的**(block structure)


解线性方程组的方法可以分为两类


1. **直接法**(direct method)
经过有限步四则运算可球的方程组准确解的方法
2. **迭代法**(iterative method)
从一个近似值出发，构造某种算法，使其逐步接近准确解


大多科学计算应用经过建模和数值离散后，都可归结为如下两种形式方程组的求解：
**方程组形式**


{a11x1\+a12x2\+⋯\+a1nxn\=b1,a21x1\+a22x2\+⋯\+a2nxn\=b2,⋯an1x1\+an2x2\+⋯\+annxn\=bn**矩阵形式**


\[a11a12⋯a1na21a22⋯a2n⋯an1an2⋯ann]\[x1x2⋮xn]\=\[b1b2⋮bn]
> Ax\=b有唯一解⟺A非奇异


### C\+\+中的线性方程组




---


在线性代数中，一矩阵的尺寸通常称为**阶数**(order)或**维度**(dimension)。以下示例代码在主函数中定义了稀疏矩阵A，常向量b和解向量x。


在**`Eigen`**库中，可以采用`Eigen::MatrixXd`表示矩阵类型，采用`Eigen::VectorXd`表示向量类型。矩阵和向量的尺寸可以在创建时进行设定。



> 需要注意的是，**`Eigen`**库中`Eigen::VectorXd`默认为列向量，如果需要将其作为行向量进行运算，需要在使用时进行转置，例如：`X.transpose()`



> 即使没有硬性的要求，但还是建议读者使用`const size_t`类型的变量单独存储矩阵的尺寸，这将使得代码维护变得更容易。



```


|  | #include |
| --- | --- |
|  | #include |
|  |  |
|  | int main() { |
|  | // 矩阵的阶数 |
|  | const size_t order = 6; |
|  |  |
|  | // 定义系数矩阵 A |
|  | Eigen::MatrixXd A(order, order); // 指定尺寸为 order * order |
|  | // 定义常向量 b |
|  | Eigen::VectorXd b(order);  // 指定尺寸为 order * 1 |
|  | // 定义解向量 x |
|  | Eigen::VectorXd x(order); // 指定尺寸为 order * 1 |
|  | } |


```

采用直接法求解线性方程组的求解器通常包含三个输入，即：系数矩阵A、常向量b和解向量x。
在进行求解前，首先应当检查输入是否符合求解器要求，例如针对上三角矩阵的求解器需要检查系数矩阵是否为上三角矩阵；一般的，输入应满足三个要求：


1. 系数矩阵A为方阵
2. 系数矩阵A的行数等于常向量b的行数
3. 系数矩阵A的列数等于解向量x的行数



> 矩阵的行数可以通过`.rows()`方法得到
> 矩阵的列数可以通过`.cols()`方法得到



> 该方法对于向量同样适用，特别的，**`Eigen`**库中向量的列数总是`1`


以下给出参考的实现：



```


|  | void size_check(const Eigen::MatrixXd& A, |
| --- | --- |
|  | const Eigen::VectorXd& b, |
|  | const Eigen::VectorXd& x) |
|  | { |
|  | // 检查A是否为方阵 |
|  | if (A.rows() != A.cols()) { |
|  | throw std::invalid_argument("Error: The coefficient matrix of the system of equations is not a square matrix."); |
|  | } |
|  | // 检查系数矩阵A的尺寸是否与常向量b的尺寸匹配 |
|  | if (A.rows() != b.rows()) { |
|  | throw std::invalid_argument("Error: The order of the coefficient matrix A does not match the order of the constant vector b."); |
|  | } |
|  | // 检查系数矩阵A的尺寸是否与解向量x的尺寸匹配 |
|  | if (A.cols() != x.rows()) { |
|  | throw std::invalid_argument("Error: The order of the coefficient matrix A does not match the order of the solution vector x."); |
|  | } |
|  | } |
|  |  |
|  | void solve(const Eigen::MatrixXd& A, |
|  | const Eigen::VectorXd& b, |
|  | Eigen::VectorXd& x) |
|  | { |
|  | // 检查尺寸是否合适 |
|  | size_check(A, b, x); |
|  |  |
|  | // 求解 |
|  | // ... |
|  | } |


```

在实际实现时有几个应注意的细节


**为什么不将解向量x作为输出？**
将解向量x作为输出的函数的使用方式为：`ans=solve(A,b)`，如果返回值的尺寸与变量`ans`的尺寸不一致则会导致程序错误。为了避免该问题，必须在创建变量`ans`时设置尺寸，并在求解前检查尺寸，伪代码如下



```


|  | Eigen::VetorXd x(order); |
| --- | --- |
|  |  |
|  | if (x.rows() == A.cols()) { // 尺寸检查 |
|  | x = solve(A,b); |
|  | } |


```

显然，形如`ans=solve(A,b,x)`的求解器更为易用，其类型检查可以在函数内部完成，这带来了更好的封装性、可维护性。


**在必要的时候添加`&`和`const`关键字**
在传递函数参数时，`&`关键字表明了该传参方式为引用传参，区别于普通传参，引用传参方式使得函数无需在其内部拷贝一个副本，而是可以直接在原变量上进行操作。无需拷贝副本显著降低了程序的性能开销。
对于普通传参，`const`关键字表明内部拷贝的副本为常变量。对于引用传参，`const`关键字表明该函数不具有修改该变量的权限，只具备读取（访问）的权限。


## 三角方程组


### 下三角方程组




---


\[a11a21a22⋮⋱an1an2⋯ann]\[x1x2⋮xn]\=\[b1b2⋮bn]解法：**前代法**(Forward substitution)


{x1\=b1/a11x2\=(b2−a21x1)/a22⋯xi\=(bi−∑j\=1i−1aijxj)/aii,i\=1,2,⋯,n**下三角矩阵判断**
**`Eigen`**库并没有提供直接的判断矩阵是否为下三角矩阵的方法，因此采用了如下的判断方法：


1. 首先提取矩阵的严格上三角部分（不包含对角线）
2. 判断其是否全部为零，如果严格上三角部分全部为零，那么其为下三角矩阵


**前代法求解**


1. 检查输入尺寸是否匹配
2. 判断系数矩阵是否为下三角矩阵
3. 采用前代法求解。


(2\.1\)xi\=(bi−∑j\=1i−1aijxj)/aii,i\=1,2,⋯,n外层循环用于遍历解向量x的每个元素，从下标`0`开始，遍历至下标`n-1`结束。循环内部分布实现式(2\.1)的计算，对于求和部分，嵌套内层循环实现。


**矩阵/向量元素访问**
在访问矩阵/向量的元素时元素，采用括号运算符进行访问。



```


|  | #include "check.h" |
| --- | --- |
|  |  |
|  | bool isLowerTriangular(const Eigen::MatrixXd& A) { |
|  | // 获取矩阵的严格上三角部分（不包括对角线） |
|  | Eigen::MatrixXd upperTriangularPart = A.triangularView(); |
|  |  |
|  | // 检查严格上三角部分是否全为零 |
|  | return upperTriangularPart.isZero(); |
|  | } |
|  |  |
|  | void forward_substitution(const Eigen::MatrixXd& A, |
|  | const Eigen::VectorXd& b, |
|  | Eigen::VectorXd& x) |
|  | { |
|  | // 检查尺寸是否匹配 |
|  | size_check(A, b, x); |
|  | // 判断系数矩阵是否为下三角矩阵 |
|  | if (!isLowerTriangular(A)) { |
|  | throw std::invalid_argument("Error: The matrix is not lower triangular."); |
|  | } |
|  |  |
|  | for (size_t i = 0; i < A.rows(); ++i) { |
|  | x(i) = b(i); |
|  | for (size_t j = 0; j + 1 <= i; ++j) { // j < i - 1 |
|  | x(i) -= A(i, j) * x(j); |
|  | } |
|  | x(i) /= A(i, i); |
|  | } |
|  | } |


```

**注意事项**



> 应当注意C\+\+中的数组索引是从`0`开始的，**`Eigen`**库也沿用了这一习惯。



> 在求和∑j\=1i−1aijxj的实现中，很容易错误的使用`j<=i-1`作为循环的终止条件，这实际上有一个风险，当`i=0`的时候，`i-1`并不是\-1，而是最大的`size_t`类型的数，这将导致终止条件错误，因此，应当用`j+1<=i`


### 上三角方程组




---


\[a11a12⋯a1na22⋯a2n⋱⋮ann]\[x1x2⋮xn]\=\[b1b2⋮bn]解法：**回代法(Back substitution)**


{xn\=bn/annxn−1\=(bn−1−an−1,nxn)/an−1,n−1⋯xi\=(bi−∑j\=1i−1aijxj)/aii,i\=n,n−1,⋯,1**上三角矩阵判断**
**`Eigen`**库并没有提供直接的判断矩阵是否为上三角矩阵的方法，因此采用了如下的判断方法：


1. 首先提取矩阵的严格下三角部分（不包含对角线）
2. 判断其是否全部为零，如果严格下三角部分全部为零，那么其为上三角矩阵


**回代法求解**


1. 检查输入尺寸是否匹配
2. 判断系数矩阵是否为上三角矩阵
3. 采用回代法求解。


(2\.2\)xi\=(bi−∑j\=1i−1aijxj)/aii,i\=n,n−1,⋯,1外层循环用于遍历解向量x的每个元素，从下标`n-1`开始，遍历至下标`0`结束。循环内部分布实现式(2\.2)的计算，对于求和部分，嵌套内层循环实现。



```


|  | bool isUpperTriangular(const Eigen::MatrixXd& A) { |
| --- | --- |
|  | // 获取矩阵的严格下三角部分（不包括对角线） |
|  | Eigen::MatrixXd lowerTriangularPart = A.triangularView(); |
|  |  |
|  | // 检查严格下三角部分是否全为零 |
|  | return lowerTriangularPart.isZero(); |
|  | } |
|  |  |
|  | void back_substitution(const Eigen::MatrixXd& A, |
|  | const Eigen::VectorXd& b, |
|  | Eigen::VectorXd& x) |
|  | { |
|  | // 检查尺寸是否匹配 |
|  | size_check(A, b, x); |
|  | // 判断系数矩阵是否为上三角矩阵 |
|  | if (!isUpperTriangular(A)) { |
|  | throw std::invalid_argument("Error: The matrix is not upper triangular."); |
|  | } |
|  |  |
|  | size_t n = A.rows(); |
|  | for (size_t i = n - 1; i != size_t(-1); --i) { // i != -1 |
|  | x(i) = b(i); |
|  | for (size_t j = i + 1; j <= n - 1; ++j) { |
|  | x(i) -= A(i, j) * x(j); |
|  | } |
|  | x(i) /= A(i, i); |
|  | } |
|  | } |


```

**注意事项**



> 外层循环的遍历是从下标`n-1`开始，遍历至下标`0`结束；一般习惯性的写法是，以`i>=0`作为截止条件，但应当注意，`size_t`类型是非负的，事实上，对于`size_t`类型的变量，当其值为`0`时再做`-1`，其值为`size_t(-1)`，因此，可以采用`i!=size_t(-1)`作为截止条件


## 高斯消元法


### 一般高斯消元法




---


**高斯消元法**（Gaussian Elimination）是一种用于求解线性方程组的经典方法。它通过逐步消去未知数，将方程组化为上三角形式，然后通过回代法求解未知数。高斯消元法主要分为两个步骤：**前向消元**和后向回代，本文中将以前向消元为例展开讨论。


**前向消元（Forward Elimination）**
前向消元法是从第一列开始，通过一些列的行变换，逐渐将原矩阵变换为一个上三角矩阵。假定矩阵的尺寸为N∗N，那么高斯消元法需要进行N−1次，在第i时执行如下操作：


1. 选择主元：选择第i列的元素Ai,i作为主元
2. 消去操作：通过将第i行的适当倍数加到其他行，使得当前列的其它元素变为零。


消去操作的公式如下：


(3\.1\){mik\=aik(k)/akk(k)aij(k\+1)\=aij(k)−m⋅akj(k)bi(k\+1)\=bi(k)−m⋅bk(k)k\=1,2,…,n−1i,j\=k\+1\.…,n矩阵的第一步消元过程可以参考以下公式：


\[a11(1)a12(1)⋯a1n(1)b1(1)a21(1)a22(1)⋯a2n(1)b2(1)⋯⋮an1(1)an2(1)⋯ann(1)bn(1)]⟶\[a11(1)a12(1)⋯a1n(1)b1(1)0a22(2)⋯a2n(2)b2(2)⋯⋮0an2(2)⋯ann(2)bn(2)]在下述程序中，采样行向量相减的方式实现高斯消元法，相较于逐个元素相减，代码更简洁易懂，易维护。



```


|  | void simple_gauss_elimination(Eigen::MatrixXd& A, Eigen::VectorXd& b) { |
| --- | --- |
|  | // 检查尺寸是否匹配 |
|  | size_check(A, b); |
|  |  |
|  | size_t n = A.rows(); |
|  | // 逐步消元为上三角矩阵 |
|  | for (size_t k = 0; k < n - 1; ++k) { |
|  | // 提取矩阵的第k行 |
|  | Eigen::VectorXd temp = A.row(k); |
|  | // 将第i列索引大于i的元素消为0 |
|  | for (size_t i = k + 1; i < n; ++i) { |
|  | // 计算比值 |
|  | double m = A(i, k) / A(k, k); |
|  | // 消元 |
|  | A.row(i) -= m * temp; |
|  | b(i) -= m * b(k); |
|  | } |
|  | } |
|  | } |


```

### 改进的高斯消元法




---


若akk(k)→0，则m\=aik(k)/akk(0)→∞，此时直接用高斯消元法求解线性方程组是会由于舍入误差的扩大，而导致解失真。


因此在原高斯消元法的基础上，可以做改进，新增主元的选择过程，该方法称为**列主元法**，具体流程如下：


1. 寻找第k列中第k行到第n行最大的元素，记为ajk


pivot\=maxk≤i≤n\|A(i,k)\|2. 将第j行与第k行交换
3. 进行高斯消元法



```


|  | void gauss_elimination(Eigen::MatrixXd& A, Eigen::VectorXd& b) { |
| --- | --- |
|  | // 检查尺寸是否匹配 |
|  | size_check(A, b); |
|  |  |
|  | size_t n = A.rows(); |
|  | // 逐步消元为上三角矩阵 |
|  | for (size_t k = 0; k < n; ++k) { |
|  | // 选择主元 |
|  | size_t j = k; |
|  | double max = abs(A(j, k)); |
|  | for (size_t i = k + 1; i < n; ++i) { |
|  | double d = abs(A(i, k)); |
|  | if (d > max) { // 选择绝对值最大的元素 |
|  | j = i; max = d; |
|  | } |
|  | } |
|  | // 交换主元 |
|  | if (j != k) { |
|  | Eigen::VectorXd temp = A.row(j); |
|  | A.row(j) = A.row(k); |
|  | A.row(k) = temp; |
|  | double temp_b = b(j); |
|  | b(j) = b(k); |
|  | b(k) = temp_b; |
|  | } |
|  | // 将第i列索引大于i的元素消为0 |
|  | for (size_t i = k + 1; i < n; ++i) { |
|  | // 计算比值 |
|  | double m = A(i, k) / A(k, k); |
|  | // 消元 |
|  | A.row(i) -= m * A.row(k); |
|  | b(i) -= m * b(k); |
|  | } |
|  | } |
|  | } |


```

**注意事项**



> 对方程Ax\=b的系数矩阵A和常向量b同时做行变换时，方程的解x不变。


### 基于高斯消元法的一般线性方程求解




---


对于一般的线性方程组，可以先用高斯消元法将系数矩阵转化为上三角矩阵，再通过回代法求解。



```


|  | void gauss_solve(Eigen::MatrixXd A, |
| --- | --- |
|  | Eigen::VectorXd b, |
|  | Eigen::VectorXd& x) |
|  | { |
|  | // 检查尺寸是否匹配 |
|  | size_check(A, b, x); |
|  | // 高斯消元法转为上三角矩阵 |
|  | gauss_elimination(A, b); |
|  | // 通过回代法求解 |
|  | back_substitution(A, b, x); |
|  | } |


```

**注意事项**



> 切忌舍本逐末，虽然添加引用修饰符可以一定程度上提升性能，但是这会导致稀疏矩阵A和常向量b被修改，而用户往往容易忽略这一点，因此为了保证**安全性**，此处**不使用引用传参**。



> 截止到目前，对系数矩阵A为下三角形矩阵的线性方程组有两种求解方法，一种是采用**前代法**，一种是采用**高斯消元结合回代法**，在附录中我们对同一组数据采用两种方法分别计算结果，进行**交叉验证**。


## 附录


### 功能测试方法




---


构建函数（方法）的测试程序流程如下：


1. 从函数（方法）的名称中提取缩写，作为名声空间的前缀
2. 定义测试函数，命名为`test()`，如果需要可以设计多个，例如：`test1()`, `test2()`
3. 实现测试函数，一般来说，有以下步骤：①生成数据，②调用方法，③打印数据以及结果
4. 在主函数中，调用该名声空间下的测试函数`test()`，一般需要使用`try-catch`结构


示例代码如下：



```


|  | namespace SMP{ |
| --- | --- |
|  | void test() { |
|  | std::cout << "Hello World!"; |
|  | } |
|  | } |
|  |  |
|  | int main() { |
|  | try{ |
|  | SMP::test(); |
|  | } |
|  | catch (const std::exception& e) { |
|  | std::cerr << "Error: " << e.what() << std::endl; |
|  | } |
|  | } |


```

在后续的附录内容中，将省略`main`函数的设计，读者只需按照上述方法调用即可。


### 前代法测试




---



```


|  | namespace FWD{ |
| --- | --- |
|  | // test for forward_substitution() |
|  | void test() { // 矩阵的阶数 |
|  | const size_t order = 5; |
|  |  |
|  | // 定义系数矩阵 A |
|  | Eigen::MatrixXd A(order, order); |
|  | // 定义常向量 b |
|  | Eigen::VectorXd b(order); |
|  | // 定义解向量 x |
|  | Eigen::VectorXd x(order); |
|  |  |
|  | // 设置矩阵为随机数 |
|  | A.setRandom(); |
|  | b.setRandom(); |
|  |  |
|  | // 处理为方便手算的数字 |
|  | A = (1.5 + A.array()) * 2;; |
|  | b *= 10; |
|  | A = A.array().round().matrix(); |
|  | b = b.array().round().matrix(); |
|  |  |
|  | // 将严格上三角部分设置为零，使其成为下三角矩阵 |
|  | A.triangularView().setZero(); |
|  |  |
|  | // 前代法 |
|  | forward_substitution(A, b, x); |
|  |  |
|  | // 输出结果 |
|  | std::cout << "A=\n" << A << "\n"; |
|  | std::cout << "b=\n" << b << "\n"; |
|  | std::cout << "x=\n" << x << "\n"; |
|  | } |
|  | } |


```

**效果展示**
程序的输出如下图所示（经过拼接），经检验，该计算结果正确（读者感兴趣的可以手算一下试试）。
![description](https://img2024.cnblogs.com/blog/3320410/202411/3320410-20241127134719109-1589879938.png)


### 回代法测试




---



```


|  | namespace BCK{ |
| --- | --- |
|  | // test for back_substitution() |
|  | void test() { // 矩阵的阶数 |
|  | const size_t order = 5; |
|  |  |
|  | // 定义系数矩阵 A |
|  | Eigen::MatrixXd A(order, order); |
|  | // 定义常向量 b |
|  | Eigen::VectorXd b(order); |
|  | // 定义解向量 x |
|  | Eigen::VectorXd x(order); |
|  |  |
|  | // 设置矩阵为随机数 |
|  | A.setRandom(); |
|  | b.setRandom(); |
|  |  |
|  | // 处理为方便手算的数字 |
|  | A = (1.5 + A.array()) * 2;; |
|  | b *= 10; |
|  | A = A.array().round().matrix(); |
|  | b = b.array().round().matrix(); |
|  |  |
|  | // 将严格下三角部分设置为零，使其成为上三角矩阵 |
|  | A.triangularView().setZero(); |
|  |  |
|  | // 回代法 |
|  | back_substitution(A, b, x); |
|  |  |
|  | // 输出结果 |
|  | std::cout << "A=\n" << A << "\n"; |
|  | std::cout << "b=\n" << b << "\n"; |
|  | std::cout << "x=\n" << x << "\n"; |
|  | } |
|  | } |


```

**效果展示**
程序的输出如下图所示（经过拼接），经检验，该计算结果正确.
![description](https://img2024.cnblogs.com/blog/3320410/202411/3320410-20241127143544259-47773379.png)


### 一般高斯消元法测试




---



```


|  | namespace S_GSE { |
| --- | --- |
|  | // test for simple_gauss_elimination |
|  | void test() { // 矩阵的阶数 |
|  | const size_t order = 5; |
|  |  |
|  | // 定义系数矩阵 A |
|  | Eigen::MatrixXd A(order, order); |
|  | // 定义常向量 b |
|  | Eigen::VectorXd b(order); |
|  |  |
|  | // 设置矩阵为随机数 |
|  | A.setRandom(); |
|  | b.setRandom(); |
|  |  |
|  | // 调整显示精度为小数点后两位 |
|  | std::cout << std::fixed << std::setprecision(2); |
|  |  |
|  | // 输出消元前矩阵 |
|  | std::cout << "A=\n" << A << "\n"; |
|  | std::cout << "b=\n" << b << "\n"; |
|  |  |
|  | // 前代法 |
|  | simple_gauss_elimination(A, b); |
|  |  |
|  | // 输出消元后矩阵 |
|  | std::cout << "A=\n" << A << "\n"; |
|  | std::cout << "b=\n" << b << "\n"; |
|  | } |
|  | } |


```

**效果展示**
程序的输出如下图所示（经过拼接），显示精度为小数点后两位；经检验，该计算结果正确.
![description](https://img2024.cnblogs.com/blog/3320410/202411/3320410-20241127161301933-1486113131.png)


### 列主元法改进的高斯消元法测试




---



```


|  | namespace GSE { |
| --- | --- |
|  | // test for simple_gauss_elimination |
|  | void test() { // 矩阵的阶数 |
|  | const size_t order = 5; |
|  |  |
|  | // 定义系数矩阵 A |
|  | Eigen::MatrixXd A(order, order); |
|  | // 定义常向量 b |
|  | Eigen::VectorXd b(order); |
|  |  |
|  | // 设置矩阵为随机数 |
|  | A.setRandom(); |
|  | b.setRandom(); |
|  |  |
|  | // 调整显示精度为小数点后两位 |
|  | std::cout << std::fixed << std::setprecision(2); |
|  |  |
|  | // 输出消元前矩阵 |
|  | std::cout << "A=\n" << A << "\n"; |
|  | std::cout << "b=\n" << b << "\n"; |
|  |  |
|  | // 前代法 |
|  | gauss_elimination(A, b); |
|  |  |
|  | // 输出消元后矩阵 |
|  | std::cout << "A=\n" << A << "\n"; |
|  | std::cout << "b=\n" << b << "\n"; |
|  | } |
|  | } |


```

程序的输出如下图所示（经过拼接），显示精度为小数点后两位；经检验，该计算结果正确.
![description](https://img2024.cnblogs.com/blog/3320410/202411/3320410-20241127194839335-1936874072.png)


### 高斯\+回代法求解




---



```


|  | namespace GS_SOLVE{ |
| --- | --- |
|  | void test1() { |
|  | const size_t order = 5; |
|  |  |
|  | Eigen::MatrixXd A(order, order); |
|  | Eigen::VectorXd b(order); |
|  | Eigen::VectorXd x(order); |
|  |  |
|  | // 设置矩阵为随机数 |
|  | A.setRandom(); |
|  | b.setRandom(); b = (1.0 + b.array()) * 5; |
|  |  |
|  | // 前代法 |
|  | gauss_solve(A, b, x); |
|  |  |
|  | // 输出结果 |
|  | std::cout << std::fixed << std::setprecision(2); |
|  | std::cout << "A=\n" << A << "\n"; |
|  | std::cout << "b=\n" << b << "\n"; |
|  | std::cout << "x=\n" << x << "\n"; |
|  | } |
|  |  |
|  | void test2() { |
|  | const size_t order = 5; |
|  |  |
|  | Eigen::MatrixXd A(order, order); |
|  | Eigen::VectorXd b(order); |
|  | Eigen::VectorXd x1(order); |
|  | Eigen::VectorXd x2(order); |
|  |  |
|  | // 设置矩阵为随机数 |
|  | A.setRandom(); |
|  | b.setRandom(); b = (1.0 + b.array()) * 5; |
|  |  |
|  | // 将上三角部分设置为零，使其成为下三角矩阵 |
|  | A.triangularView().setZero(); |
|  |  |
|  | // 高斯 |
|  | gauss_solve(A, b, x1); |
|  | // 前代法 |
|  | forward_substitution(A, b, x2); |
|  |  |
|  | // 输出结果 |
|  | std::cout << std::fixed << std::setprecision(2); |
|  | std::cout << "GS_solve:\n" << "x1=\n" << x1 << "\n"; |
|  | std::cout << "back_stt:\n" << "x2=\n" << x2 << "\n"; |
|  | } |
|  | } |


```

**测试1**
函数`GS_SOLVE::test1()`用于测试高斯求解是否能够正常工作，该程序的输出如下图所示（经过拼接），显示精度为小数点后两位；经检验，该计算结果正确.
![description](https://img2024.cnblogs.com/blog/3320410/202411/3320410-20241127202608712-196359627.png)


**测试2**
函数`GS_SOLVE::test2()`采用交叉验证法，分别采用**前代法**，一种是采用**高斯消元结合回代法**求解系数矩阵A为下三角矩阵的线性方程组，并对比计算结果；经检验，结果各方面功能正常。
![description](https://img2024.cnblogs.com/blog/3320410/202411/3320410-20241127203030853-291746863.png)


 本博客参考[FlowerCloud机场](https://hushicha.org)。转载请注明出处！
