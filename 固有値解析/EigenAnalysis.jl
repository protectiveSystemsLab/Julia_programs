#=
せん断型構造物の固有値解析サンプル
社団法人　免震構造協会　編，パッシブ制震構造　設計・施工マニュアル　第２版
４階建てテーマストラクチャー在来タイプのせん断棒モデルの解析
=#
using LinearAlgebra
using Plots
gr()                                # バックエンドの指定（デフォルトはGR）

# 単位系 t, cm , kN
g = 980.0                           # 重力加速度
n = 4                               # 自由度数4
w = [3863.1 3794.2 3731.4 4966.1]   # 各階重量(成分は1階から4階の順)
k = [1427.9 2183.1 1736.1 1406.0]   # 各階せん断剛性
h = [600 400 400 400]               # 各階階高

M  = zeros(n,n)
kd = zeros(n,n)
H  = zeros(n+1)
H[1] = 0.0
for j = 1:n
    M[j,j]  = w[j]/g                 # 質量行列
    kd[j,j] = k[j]                   # 各階せん断剛性を対角に持つ行列
    H[j+1]  = H[j]+h[j]
end

T = Matrix(I, n, n)-diagm(-1 => ones(n-1))
K = T'*kd*T

lambda = eigen(M\K)
omega = sqrt.(lambda.values)
index = sortperm(lambda.values)
# omega = omega[index]
P = 2π./omega

modeStr = ["1st" "2nd" "3rd" "4th"]
println("固有周期")
for j=1:n
    period = P[j]
    println(modeStr[j], " mode: $period")
end

U = lambda.vectors[:,index]

r = ones(n)                         # 影響係数ベクトル

Me = U'*M*U                         # 一般化質量の計算
me = diag(Me)                       # 対角成分の抽出
nu = (U'*M*r)./me                   # 刺激係数

Nu = zeros(n,n)
for j=1:n
    Nu[j,j] = nu[j]
end

nuU0 = U * Nu                        # 刺激関数

nuU = zeros(n+1,n)
nuU[2:5,1:4]=nuU0

plot(nuU,H,
    title="Mode shape",
    size=(400,400),
    label = ["1st mode" "2nd mode" "3rd mode" "4th mode"],
    xlims = (-1.5,1.5),
    ylims = (0,Inf),
    lw = 2)
