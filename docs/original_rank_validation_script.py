exec(open('/mnt/data/run_sample_size_calibration_fast.py').read().split("rows=[]\ncond=0")[0])
# Functions only via split before execution loops
cases=[(176,.10),(222,.10),(230,.10),(278,.10),(60,.30),(72,.30),(77,.30),(35,.30)]
for n,R2 in cases:
    sels=[]
    for j in range(500):
        seed=50000000+j
        rng=np.random.default_rng(seed)
        X=rng.normal(size=(n,24)); B=make_B(rng,24,24,10,R2,'equal'); Y=X@B+rng.normal(size=(n,24))
        sel,_=cv_select(X,Y,seed); sels.append(sel)
    a=np.array(sels)
    print(n,R2,'median',np.median(a),'mean',a.mean(),'q25q75',np.quantile(a,[.25,.75]),'exact',np.mean(a==10),'zero',np.mean(a==0))
