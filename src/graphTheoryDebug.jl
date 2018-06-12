include("auxiliar.jl")
module GraphTheory
    using Auxiliar
    function WalkTail(latt,pos,neigLatt;printLog=false)
        l=[pos]
        nm=deepcopy(neigLatt)
        sameNeigs=[x for x in nm[pos] if latt[x]==latt[pos] && in(pos,nm[x])]
        newPos=sameNeigs[1]
        #self loop exception
        if length(Set(sameNeigs))==1 &&  newPos==pos
            return []
        end
        push!(l,newPos)
        #deleting old links
        index=Auxiliar.Index(nm[pos],newPos)
        deleteat!(nm[pos],index)
        index=Auxiliar.Index(nm[newPos],pos)
        deleteat!(nm[newPos],index)
        newNeigs=[x for x in nm[newPos] if latt[x]==latt[newPos]]
        while ~(isempty(newNeigs) || length(newNeigs)>=2)
            pos=newPos
            sameNeigs=[x for x in nm[pos] if latt[x]==latt[pos] && in(pos,nm[x])]
            newPos=sameNeigs[1]
            push!(l,newPos)
            index=Auxiliar.Index(nm[pos],newPos)
            deleteat!(nm[pos],index)
            index=Auxiliar.Index(nm[newPos],pos)
            deleteat!(nm[newPos],index)
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
        sameNeigs=[x for x in neigLatt[pos] if latt[x]==latt[pos] && in(pos,neigLatt[x])]
        #special case of single tail
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
        m=deepcopy(latt)
        while sum(vis) < prod(size(latt))
            pos=Auxiliar.RandomPosition(latt)
            l=Tail(m,pos,neigLatt,printLog=printLog)
            if printLog
                println(pos)
                println(l)
            end
            if isempty(l)
                vis[pos]=1
            else
                for p in l
                    m[p]=2
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
        if printLog
            println()
            println(" entering simple walk")
            println()
        end
        l=[]
        initPos=pos
        nm=deepcopy(neigLatt)
        sameNeigs=[x for x in nm[pos] if latt[x]==latt[pos] && in(pos,nm[x])]
        if ~isempty(sameNeigs)
            push!(l,pos)
            newPos=rand(sameNeigs)
            #self loop exception
            if length(Set(sameNeigs))==1 &&  newPos==pos
                return [pos]
            end
            index=Auxiliar.Index(nm[pos],newPos)
            deleteat!(nm[pos],index)
            index=Auxiliar.Index(nm[newPos],pos)
            try deleteat!(nm[newPos],index)
            catch LoadError
            end
            while initPos!=newPos && ~in(newPos,l)
                pos=newPos
                sameNeigs=[x for x in nm[pos] if latt[x]==latt[pos] && in(pos,nm[x])]
                if ~isempty(sameNeigs)
                    push!(l,pos)
                    newPos=rand(sameNeigs)
                    index=Auxiliar.Index(nm[pos],newPos)
                    deleteat!(nm[pos],index)
                    index=Auxiliar.Index(nm[newPos],pos)
                    try deleteat!(nm[newPos],index)
                    catch LoadError
                        continue
                    end
                else
                    push!(l,pos)
                    break
                end
            end
            if (initPos==newPos)
                return l
            else
                index=Auxiliar.Index(l,newPos)
                return l[index:end]
            end
        else
            return [pos]
        end
        if printLog
            println()
            println(" leaving simple walk")
            println()
        end
    end
    function WalkComplicatedPath(latt,pos,neigLatt;printLog=false) 
        if printLog
            println()
            println(" entering complicated walk")
            println()
        end
        l=WalkSimplePath(latt,pos,neigLatt;printLog=printLog)
        nm=deepcopy(neigLatt)
        if length(l)==1
            l=[x for x in l]
            nm[pos]=[]
            return ([l[1],l[1],l[1]],nm)
        end
        if ~isempty(l)
            b=[] 
            for p in l
                #println("infinite loop vol2")
                n=[neig for neig in nm[p] if ~in(neig,l) && in(p,nm[neig])]
                x=[neig for neig in n if (latt[p]==latt[neig])]        
                for t in n
                    index=Auxiliar.Index(nm[t],p)
                    deleteat!(nm[t],index)
                end
                append!(b,x)
                nm[p]=[]
            end
            #=
            Exception for single points between cycles
            
                Example (4,3) in
            [
                0 0 0 0 0 ;
                0 1 1 1 0 ;
                0 1 0 1 0 ;
                0 1 1 1 0 ; 
                0 1 0 1 0 ;
                0 1 1 1 0 ;
            ]
            =#
            if length(Set(b))==1 && length(b)>1
                p=b[1]
                push!(nm[p],p)
            else
                b=Set(b)
                if length(b)==2
                    b=[x for x in b]
                    append!(nm[b[1]],[b[2],b[2]])
                    append!(nm[b[2]],[b[1],b[1]])
                else
                    for p in b
                        #println("infinite loop vol3")
                        k=[a for a in b if a!=p]
                        append!(nm[p],k)
                        nm[p]=[x for x in Set(nm[p])]
                    end
                end
            end
            return (Set(l),nm)
        else
            return ([pos],nm)
        end
        if printLog
            println()
            println(" leaving complicated walk")
            println()
        end
    end


    function RemoveElements!(latt,l;printLog=false)
        for pos in l
            latt[pos]=2
        end
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
    function RemoveDuplicates!(neigLatt)
        for pos in CartesianRange(size(neigLatt))
            neigLatt[pos]=[a for a in Set(neigLatt[pos])]
        end
    end

    function SearchCycles(latt,neigLatt;printLog=false)
        if printLog
            println()
            println(" entering cycle search")
            println()
        end
        m=copy(latt)
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
                x=WalkComplicatedPath(m,pos,nm,printLog=printLog)
                l=x[1]
                nm=x[2]
                RemoveDuplicates!(nm)
                RemoveElements!(m,l)
                m3=deepcopy(m)
                m=RemoveTails(m,nm)
                for p in l
                    vis[p]=1
                end
                if length(l)>=2
                    if printLog
                        println("cycle")
                        println(l)
                        println("old matrix before tailing")
                        println(m3)
                        println("new matrix")
                        println(m)
                        println("neigbors")
                        println(nm)
                    end
                    append!(fin,l)
                end
            end
        end
        if printLog
            println()
            println(" leaving cycle search")
            println()
        end
        return (Set(fin),x[2])
    end

    function Cycles(latt,neigLatt;printLog=false)
        d=RemoveTails(RemoveOthers(deepcopy(latt),-1),neigLatt)
        l=SearchCycles(d,neigLatt;printLog=printLog)
        return l[1]
    end
    a=Array{Int8,2}([-1 1 1 1 -1 1 -1 -1; -1 1 1 1 -1 1 -1 1; 1 1 1 1 1 1 1 -1; 1 -1 -1 1 1 -1 -1 -1; 1 -1 1 -1 1 -1 1 -1; -1 -1 1 1 -1 -1 -1 1; -1 -1 1 1 -1 -1 -1 1; 1 -1 -1 1 1 -1 1 -1])
    neigLatt=Auxiliar.NeighborIndexLattice(a,Auxiliar.SquareLatticeNeighborsIndex)
    c=Cycles(a,neigLatt,printLog=false)
    println("a")
    println(a)
    println("c")
    println(c)
    for i in 1:200
        @time b=Cycles(a,neigLatt,printLog=true)
        if c!=b
            println(b)
            println("fail")
        end
    end
end