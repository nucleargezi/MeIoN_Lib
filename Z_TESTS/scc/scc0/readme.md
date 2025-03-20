### scc test1

给定 $n$ 个点 $m$ 条边的一张有向图，求它的scc

#### input 
$n$ $m$  
$x_1$ $y_1$  
$x_2$ $y_2$  
...  
$x_m$ $y_m$

#### output

输出的答案是一个异或和  

设图中有 $k$ 个强连通分量（SCC），记为 $S_1, S_2, \ldots, S_k$。对于每个强连通分量 $S_i$：

- $|S_i|$ 表示该 SCC 的大小（包含的节点数）
- $\bigoplus\limits_{v \in S_i} v$ 表示该 SCC 内所有节点编号的异或和

则最终的表达式可以表示为：

$$
k \oplus \left( \bigoplus_{i=1}^{k} \left( |S_i| \oplus \left( \bigoplus_{v \in S_i} v \right) \right) \right)
$$

其中：
- $\oplus$ 表示异或运算
- 外层第一个 $k$ 是 SCC 的总数
- 内层嵌套的异或运算先计算每个 SCC 的「大小 异或 节点编号异或和」
- 最后将所有 SCC 的计算结果进行异或，再与 SCC 总数 $k$ 异或