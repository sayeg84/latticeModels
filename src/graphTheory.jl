include("auxiliar.jl")
module GraphTheory
    using Auxiliar
    function Modl(a::Int64,b::Int64)
        x=mod(a,b)
        if x==0
            return b
        else
            return x
        end
    end
    function MtoV(m;printLog=false)
        v=[]
        for i in 1:size(m)[1]
            for j in 1:size(m)[2]
                push!(v,m[i,j])
            end
        end
        return v
    end
    function VtoM(v,text=false)
        s=1
        try s=convert(Int64,sqrt(length(v))) 
        catch InexactError; error("Vector is not perfect square")
        end
        m=zeros(s,s)
        for i in 1:s
            for j in 1:s 
                m[i,j]=v[j+(i-1)*s]
            end
        end
        return m
    end
    function MtoVpos(p,size)
        return p[2]+size[2]*(p[1]-1)
    end
    function VtoMpos(p,size)
        i=ceil(p/size[2])
        i=convert(Int64,i)
        j=Modl(p,size[2])
        return [i,j]
    end
    function ConstructAdjacencyMatrix(m;printLog=false)
        x=size(m)[1]
        y=size(m)[2]
        fin=zeros(x*y,x*y)
        for i in 1:x 
            for j in 1:y 
                if m[i,j]==m[Modl(i+1,x),j]
                    fin[MtoVpos([i,j],[x,y]),MtoVpos([Modl(i+1,x),j],[x,y])]=1
                end
                if m[i,j]==m[Modl(i-1,x),j]
                    fin[MtoVpos([i,j],[x,y]),MtoVpos([Modl(i-1,x),j],[x,y])]=1
                end
                if m[i,j]==m[i,Modl(j+1,y)]
                    fin[MtoVpos([i,j],[x,y]),MtoVpos([i,Modl(j+1,y)],[x,y])]=1
                end
                if m[i,j]==m[i,Modl(j-1,y)]
                    fin[MtoVpos([i,j],[x,y]),MtoVpos([i,Modl(j-1,y)],[x,y])]=1
                end
            end
        end
        return fin
    end
    function RemoveOthers(latt,oth;printLog=false)
        m=deepcopy(latt)
        for pos in CartesianRange(size(latt))
            if m[pos]==oth 
                m[pos]=2
            end
        end
        return m
    end
    function WalkTail(latt,pos,neigLatt;printLog=false)
        l=[pos]
        nm=deepcopy(neigLatt)
        sameNeigs=[x for x in nm[pos] if latt[x]==latt[pos]]
        newPos=sameNeigs[1]
        push!(l,newPos)
        i=Auxiliar.Index(nm[newPos],pos)
        deleteat!(nm[newPos],i)
        newNeigs=[x for x in nm[newPos] if latt[x]==latt[newPos]]
        while ~(isempty(newNeigs) || length(newNeigs)>=2)
            pos=newPos
            sameNeigs=[x for x in nm[pos] if latt[x]==latt[pos]]
            newPos=sameNeigs[1]
            push!(l,newPos)
            i=Auxiliar.Index(nm[newPos],pos)
            deleteat!(nm[newPos],i)
            newNeigs=[x for x in nm[newPos] if latt[x]==latt[newPos]]
        end
        if length(l)>=2
            deleteat!(l,length(l))
        end
        return l
    end
    function Tail(latt,pos,neigLatt;printLog=false)
        if (latt[pos]==-1 || latt[pos] == 0 || latt[pos] == 2 )
            return []
        end
        sameNeigs=[x for x in neigLatt[pos] if latt[x]==latt[pos]]
        if isempty(sameNeigs)
            return [pos]
        elseif length(sameNeigs)>1
            return []
        else
            return WalkTail(latt,pos,neigLatt,printLog=printLog)
        end
    end

    function RemoveTails(latt,neigLatt;printLog=false)
        vis=zeros(size(latt))
        fin=deepcopy(latt)
        while sum(vis) < prod(size(latt))
        pos=Auxiliar.RandomPosition(latt)
            l=Tail(latt,pos,neigLatt,printLog=printLog)
            if isempty(l)
                vis[pos]=1
            else
                for p in l
                    fin[p]=2
                    vis[p]=1
                end
            end
        end
        return fin
    end



    function WalkSimplePath(latt,pos,neigLatt;printLog=false)
        if (latt[pos]==0 || latt[pos]==-1 || latt[pos]==2)
            return []
        end
        l=[]
        initPos=pos
        if printLog
            println()
            println("enter simple walk")
            println()
        end
        nm=deepcopy(neigLatt)
        neig=nm[pos]
        sameNeigs=[x for x in neig if latt[x]==latt[pos]]
        if ~isempty(sameNeigs)
            push!(l,pos)
            newPos=rand(sameNeigs)
            index=Auxiliar.Index(nm[newPos],pos)
            deleteat!(nm[newPos],index)
            while initPos!=newPos && ~in(newPos,l)
                #println("aldo")
                #println(newPos)
                pos=newPos
                #println("linea 1")
                neig=nm[pos]
                #println("linea 2")
                sameNeigs=[x for x in neig if latt[x]==latt[pos]]
                #println("linea 3")
                if ~isempty(sameNeigs)
                    #println("linea 4")
                    push!(l,pos)
                    #println("linea 5")
                    newPos=rand(sameNeigs)
                    #println("linea 6")
                    index=Auxiliar.Index(nm[newPos],pos)
                    #println("linea 7")
                    deleteat!(nm[newPos],index)
                    #println("linea 8" )
                end
                #println("linea 9")
                #println(initPos)
                #println(newPos)
                #println(l)
                #println(in(newPos,l))
            end
            #println("sali")
            return l
        else
            return []
        end
        if printLog
            println()
            println("exit simple walk")
            println()
        end
    end
    function WalkComplicatedPath(latt,pos,neigLatt;printLog=false) 
        if printLog
            println()
            println("enter complicated walk")
            println()
        end
        l=WalkSimplePath(latt,pos,neigLatt;printLog=printLog)
        #println("salio")
        nm=deepcopy(neigLatt)
        if length(l)==1
            nm[pos]=[]
            return ([l[1],l[1],l[1]],nm)
        end
        if ~isempty(l)
            b=[]
            #println("llegó")
            for p in l
                n=[neig for neig in neigLatt[p] if ~in(neig,l)]
                x=[neig for neig in n if (latt[p]==latt[neig])]                
                for t in n
                    i=Auxiliar.Index(nm[t],p)
                    deleteat!(nm[t],i)
                    k=[a for a in n if a!=t]
                    append!(nm[t],k)
                end
                append!(b,x)
                nm[p]=[]
            end
            for p in b
                k=[a for a in b if a!=p]
                append!(nm[p],k)
                if length(b)==2 && b[1]==b[2]
                #    println("single")
                #    println()
                #    println(pos)
                #    println()   
                #    println("cycle")
                #    println(l)
                #    println()
                #
                #    println(p)
                #    println()
                #    println(nm[p])
                #    println()
                    push!(nm[p],p)
                #    println(nm[p])
                #    println()
                #    println(k)
                #    println()
                end
            end
            return (Set(l),nm)
        else
            return ([pos],nm)
        end
        if printLog
            println()
            println("exit complicated walk")
            println()
        end
    end
    function SearchCycles(latt,neigLatt;printLog=false)
        if printLog
            println()
            println("enter cycle search")
            println()
        end
        s=size(latt)
        vis=zeros(s)
        nm=deepcopy(neigLatt)
        fin=[]
        x=[]
        while sum(vis) < prod(s)
            pos=Auxiliar.RandomPosition(latt)
            if printLog
                println()
                println("visiting pos")
                println(pos)
                println()
            end
            #println(pos)
            if vis[pos]!=1
                x=WalkComplicatedPath(latt,pos,nm,printLog=printLog)
                l=x[1]
                nm=x[2]
                for p in l
                    vis[p]=1
                end
                if length(l)>=2
                    append!(fin,l)
                end
            end
        end
        if printLog
            println()
            println("exit cycle search")
            println()
        end
        return (Set(fin),x[2])
    end

    function Cycles(latt,neigLatt;printLog=false)
        d=RemoveTails(RemoveOthers(deepcopy(latt),-1),neigLatt)
        l=SearchCycles(d,neigLatt;printLog=printLog)
        return l[1]
    end
#=
    m1=[
    0 0 0 0 1 1 1 1 1 ;
    0 0 0 0 1 0 0 0 1 ;
    0 1 1 1 1 0 1 0 1 ;
    0 0 0 0 0 0 0 0 1 ;
    0 0 0 0 0 0 0 0 0 ;
    ]
    m2=[
    0 0 0 0 1 1 1 0 0 ;
    0 1 1 1 1 0 1 0 0 ;
    0 0 0 0 1 1 1 1 1 ;
    0 0 0 0 0 0 0 0 1 ;
    0 0 0 0 0 0 0 0 0 ;
    ]
    m3=[
    0 0 0 0 1 1 1 0 0 ;
    0 2 2 2 1 0 1 0 0 ;
    0 0 0 0 2 1 1 2 2 ;
    0 0 0 0 0 0 0 0 2 ;
    0 0 0 0 0 0 0 0 0 ;
    ]
    
    m4=[
    0 0 0 0 1 1 1 0 0 ;
    0 2 2 2 1 0 1 0 0 ;
    0 0 0 0 1 1 1 2 2 ;
    0 0 0 0 1 0 1 0 2 ;
    0 0 0 0 1 1 1 0 0 ;
    ]
 
    #l=RemoveTails(m2)
    #x=WalkPath([3,5],[2,5],m3)
    m3=[
    0 0 0 0 1 1 1 1 0 0 ;
    0 2 2 2 1 0 0 1 0 0 ;
    0 0 0 0 1 1 1 1 2 2 ;
    0 0 0 0 1 0 0 1 0 2 ;
    0 0 0 0 1 1 1 1 0 0 ;
    0 0 0 0 0 0 0 0 0 0 ;
    ]
    
    m3=[
    0 0 0 0 1 1 1 1 0 0 ;
    0 1 1 1 1 0 0 1 0 0 ;
    0 0 0 0 1 1 1 1 1 1 ;
    0 0 0 0 1 0 0 1 0 1 ;
    0 0 0 0 1 1 1 1 0 1 ;
    0 0 0 0 0 0 0 0 0 0 ;
    ]
    m3=[
        0 1 1 1 1 1 1 1 1 0 0 ;
        0 0 0 0 1 0 0 0 1 0 0 ;
        0 0 1 1 1 1 1 1 1 1 1 ;
        0 0 1 0 1 0 1 0 1 0 1 ;
        0 0 1 0 1 1 1 1 1 0 1 ;
        0 0 1 0 0 0 1 0 0 0 0 ;
        0 0 1 1 1 1 1 0 0 0 0 ;
        0 0 0 0 0 0 0 0 0 0 0 ;
        ]
   

       neigLatt=Auxiliar.NeighborIndexLattice(m3,Auxiliar.SquareLatticeNeighborsIndex)
       @time a=Cycles(m3,neigLatt)
       println(a)
    for i in 1:2000
        @time x=Cycles(m3,neigLatt)
        if x==a println(true) end
    end
    =#

    #@time x=WalkComplicatedPath(m3,pos,neigLatt)
    #println(x[1])
    #println(x[2][CartesianIndex((3,6))])
    #@time y=WalkComplicatedPath(m3,pos,neigLatt)
    #println(y[1])
    #println(y[2][CartesianIndex((3,6))])
    #@time z=WalkComplicatedPath(m3,pos,neigLatt)
    #println(z[1])
    #println(z[2][CartesianIndex((3,6))])
    #@time x=WalkSimplePath(m3,pos,neigLatt)
    #println(x)
    #pos=CartesianIndex((3,5))
    #@time x=WalkPath(m3,pos,neigLatt)
    #println(x)
    #=
    g=Erdos.Graph(ConstructAdjacencyMatrix(m))
    x=Erdos.conected_components(g)
    =#
end