include("auxiliar.jl")
module GraphTheory
    using Erdos,Auxiliar
    function Modl(a::Int64,b::Int64)
        x=mod(a,b)
        if x==0
            return b
        else
            return x
        end
    end
    function MtoV(m;test=false)
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
    function ConstructAdjacencyMatrix(m;test=false)
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

    function Tail(pos,mat;test=false)
        m=copy(mat)
        l=[]
        x=size(m)[1]
        y=size(m)[2]
        i=pos[1]
        j=pos[2]
        if (m[i,j]!=m[Modl(i+1,x),j] && m[i,j]!=m[Modl(i-1,x),j] && m[i,j]!=m[i,Modl(j+1,y)] && m[i,j]!=m[i,Modl(j-1,y)])
            push!(l,pos)
            return l
        elseif (m[i,j]!=m[Modl(i+1,x),j] && m[i,j]!=m[Modl(i-1,x),j] && m[i,j]!=m[i,Modl(j+1,y)])
            push!(l,pos)
            m[i,j]=2
            j=Modl(j-1,y)
            append!(l,Tail([i,j],m))
            return l
        elseif (m[i,j]!=m[Modl(i+1,x),j] && m[i,j]!=m[Modl(i-1,x),j] && m[i,j]!=m[i,Modl(j-1,y)]) 
            push!(l,pos)
            m[i,j]=2
            j=Modl(j+1,y)
            append!(l,Tail([i,j],m))
            return l
        elseif (m[i,j]!=m[Modl(i+1,x),j] && m[i,j]!=m[i,Modl(j-1,y)] && m[i,j]!=m[i,Modl(j+1,y)])
            push!(l,pos)
            m[i,j]=2
            i=Modl(i-1,x)
            append!(l,Tail([i,j],m))
            return l
        elseif (m[i,j]!=m[Modl(i-1,x),j] && m[i,j]!=m[i,Modl(j-1,y)] && m[i,j]!=m[i,Modl(j+1,y)])
            push!(l,pos)
            m[i,j]=2
            i=Modl(i+1,x)
            append!(l,Tail([i,j],m))
            return l
        else
            return l 
        end
    end
    function RemoveTails(m;test=false)
        x=size(m)[1]
        y=size(m)[2]
        vis=zeros(x,y)
        fin=copy(m)
        while sum(vis) < x*y
            i=rand(1:x)
            j=rand(1:y)
            if test
                print("rand pos: ")
                println([i,j])
            end
            if vis[i,j]==0
                l=Tail([i,j],m)
                if isempty(l)
                    vis[i,j]=1
                else
                    for pos in l 
                        if test
                            print("pos: ")
                            println(pos)
                            print("tail: ")
                            println(l)
                            println("vis")
                            println(vis)
                            println("clean")
                            println(fin)
                        end
                        vis[pos[1],pos[2]]=1
                        fin[pos[1],pos[2]]=2
                    end
                end
            end
        end
        return fin
    end
    function WalkPath(posInit,pos,mat;test=false)
        m=copy(mat)
        x=size(m)[1]
        y=size(m)[2]
        i=pos[1]
        j=pos[2]
        l=[[i,j]]
        if (i==posInit[1] && j==posInit[2])
            return l
        end
        if (m[i,j]==m[Modl(i+1,x),j] && m[i,j]==m[Modl(i-1,x),j] && m[i,j]!=m[i,Modl(j+1,y)] && m[i,j]!=m[i,Modl(j-1,y)])
            r=rand()
            if r<0.5
                m[posInit[1],posInit[2]]=1
                m[i,j]=2
                i=Modl(i+1,x)
                append!(l,WalkPath(posInit,[i,j],m))
            else
                m[posInit[1],posInit[2]]=1
                m[i,j]=2
                i=Modl(i-1,x)
                append!(l,WalkPath(posInit,[i,j],m))
            end
        elseif (m[i,j]==m[Modl(i+1,x),j] && m[i,j]!=m[Modl(i-1,x),j] && m[i,j]==m[i,Modl(j+1,y)] && m[i,j]!=m[i,Modl(j-1,y)])
            r=rand()
            if r<0.5
                m[posInit[1],posInit[2]]=1
                m[i,j]=2
                i=Modl(i+1,x)
                append!(l,WalkPath(posInit,[i,j],m))
            else
                m[posInit[1],posInit[2]]=1
                m[i,j]=2
                j=Modl(j+1,y)
                append!(l,WalkPath(posInit,[i,j],m))
            end
        elseif (m[i,j]==m[Modl(i+1,x),j] && m[i,j]!=m[Modl(i-1,x),j] && m[i,j]!=m[i,Modl(j+1,y)] && m[i,j]==m[i,Modl(j-1,y)])
            r=rand()
            if r<0.5
                m[posInit[1],posInit[2]]=1
                m[i,j]=2
                i=Modl(i+1,x)
                append!(l,WalkPath(posInit,[i,j],m))
            else
                m[posInit[1],posInit[2]]=1
                m[i,j]=2
                j=Modl(j-1,y)
                append!(l,WalkPath(posInit,[i,j],m))
            end
        elseif (m[i,j]!=m[Modl(i+1,x),j] && m[i,j]==m[Modl(i-1,x),j] && m[i,j]==m[i,Modl(j+1,y)] && m[i,j]!=m[i,Modl(j-1,y)])
            r=rand()
            if r<0.5
                m[posInit[1],posInit[2]]=1
                m[i,j]=2
                i=Modl(i-1,x)
                append!(l,WalkPath(posInit,[i,j],m))
            else
                m[posInit[1],posInit[2]]=1
                m[i,j]=2
                j=Modl(j+1,y)
                append!(l,WalkPath(posInit,[i,j],m))
            end
        elseif (m[i,j]!=m[Modl(i+1,x),j] && m[i,j]==m[Modl(i-1,x),j] && m[i,j]!=m[i,Modl(j+1,y)] && m[i,j]==m[i,Modl(j-1,y)])
            r=rand()
            if r<0.5
                m[posInit[1],posInit[2]]=1
                m[i,j]=2
                i=Modl(i-1,x)
                append!(l,WalkPath(posInit,[i,j],m))
            else
                m[posInit[1],posInit[2]]=1
                m[i,j]=2
                j=Modl(j-1,y)
                append!(l,WalkPath(posInit,[i,j],m))
            end
        elseif (m[i,j]!=m[Modl(i+1,x),j] && m[i,j]!=m[Modl(i-1,x),j] && m[i,j]==m[i,Modl(j+1,y)] && m[i,j]==m[i,Modl(j-1,y)])
            r=rand()
            if r<0.5
                m[posInit[1],posInit[2]]=1
                m[i,j]=2
                j=Modl(j+1,y)
                append!(l,WalkPath(posInit,[i,j],m))
            else
                m[posInit[1],posInit[2]]=1
                m[i,j]=2
                j=Modl(j-1,y)
                append!(l,WalkPath(posInit,[i,j],m))
            end  
        elseif (m[i,j]==m[Modl(i+1,x),j] && m[i,j]!=m[Modl(i-1,x),j] && m[i,j]!=m[i,Modl(j+1,y)] && m[i,j]!=m[i,Modl(j-1,y)])
            m[posInit[1],posInit[2]]=1
            m[i,j]=2
            i=Modl(i+1,x)
            append!(l,WalkPath(posInit,[i,j],m))
        elseif (m[i,j]!=m[Modl(i+1,x),j] && m[i,j]==m[Modl(i-1,x),j] && m[i,j]!=m[i,Modl(j+1,y)] && m[i,j]!=m[i,Modl(j-1,y)])
            m[posInit[1],posInit[2]]=1
            m[i,j]=2
            i=Modl(i-1,x)
            append!(l,WalkPath(posInit,[i,j],m))
        elseif (m[i,j]!=m[Modl(i+1,x),j] && m[i,j]!=m[Modl(i-1,x),j] && m[i,j]==m[i,Modl(j+1,y)] && m[i,j]!=m[i,Modl(j-1,y)])
            m[posInit[1],posInit[2]]=1
            m[i,j]=2
            j=Modl(j+1,y)
            append!(l,WalkPath(posInit,[i,j],m))
        elseif (m[i,j]!=m[Modl(i+1,x),j] && m[i,j]!=m[Modl(i-1,x),j] && m[i,j]!=m[i,Modl(j+1,y)] && m[i,j]==m[i,Modl(j-1,y)])
            m[posInit[1],posInit[2]]=1
            m[i,j]=2
            j=Modl(j-1,y)
            append!(l,WalkPath(posInit,[i,j],m))
        end
        return l
    end
    function CheckCycle(posInit,mat;test=false)
        m=copy(mat)
        x=size(m)[1]
        y=size(m)[2]
        i=posInit[1]
        j=posInit[2]
        if (m[i,j]==m[Modl(i+1,x),j] && m[i,j]==m[Modl(i-1,x),j] && m[i,j]==m[i,Modl(j+1,y)] && m[i,j]==m[i,Modl(j-1,y)])
            m[i,j]=2
            r=rand()
            if r<0.25
                l=WalkPath([i,j],[i,Modl(j+1,y)],m)
            elseif r<0.5
                l=WalkPath([i,j],[i,Modl(j-1,y)],m)
            elseif r<0.75
                l=WalkPath([i,j],[Modl(i+1,x),j],m)
            else
                l=WalkPath([i,j],[Modl(i-1,x),j],m)
            end
            if (posInit==l[end])
                return l
            else
                return []
            end
        elseif (m[i,j]!=m[Modl(i+1,x),j] && m[i,j]==m[Modl(i-1,x),j] && m[i,j]==m[i,Modl(j+1,y)] && m[i,j]==m[i,Modl(j-1,y)])
            m[i,j]=2
            r=rand()
            if r<1/3
                l=WalkPath([i,j],[Modl(i-1,x),j],m)
            elseif r<2/3
                l=WalkPath([i,j],[i,Modl(j+1,y)],m)
            else
                l=WalkPath([i,j],[i,Modl(j-1,y)],m)
            end
            if (posInit==l[end])
                return l
            else
                return []
            end
        elseif (m[i,j]==m[Modl(i+1,x),j] && m[i,j]!=m[Modl(i-1,x),j] && m[i,j]==m[i,Modl(j+1,y)] && m[i,j]==m[i,Modl(j-1,y)])
            m[i,j]=2
            r=rand()
            if r<1/3
                l=WalkPath([i,j],[Modl(i+1,x),j],m)
            elseif r<2/3
                l=WalkPath([i,j],[i,Modl(j+1,y)],m)
            else
                l=WalkPath([i,j],[i,Modl(j-1,y)],m)
            end
            if (posInit==l[end])
                return l
            else
                return []
            end
        elseif (m[i,j]==m[Modl(i+1,x),j] && m[i,j]==m[Modl(i-1,x),j] && m[i,j]!=m[i,Modl(j+1,y)] && m[i,j]==m[i,Modl(j-1,y)])
            m[i,j]=2
            r=rand()
            if r<1/3
                l=WalkPath([i,j],[Modl(i+1,x),j],m)
            elseif r<2/3
                l=WalkPath([i,j],[Modl(i-1,x),j],m)
            else
                l=WalkPath([i,j],[i,Modl(j-1,y)],m)
            end
            if (posInit==l[end])
                return l
            else
                return []
            end
        elseif (m[i,j]==m[Modl(i+1,x),j] && m[i,j]==m[Modl(i-1,x),j] && m[i,j]==m[i,Modl(j+1,y)] && m[i,j]!=m[i,Modl(j-1,y)])
            m[i,j]=2
            r=rand()
            if r<1/3
                l=WalkPath([i,j],[Modl(i+1,x),j],m)
            elseif r<2/3
                l=WalkPath([i,j],[Modl(i-1,x),j],m)
            else
                l=WalkPath([i,j],[i,Modl(j+1,y)],m)
            end
            if (posInit==l[end])
                return l
            else
                return []
            end
        elseif (m[i,j]==m[Modl(i+1,x),j] && m[i,j]==m[Modl(i-1,x),j])
            m[i,j]=2
            r=rand()
            if r<0.5
                l=WalkPath([i,j],[Modl(i+1,x),j],m)
            else
                l=WalkPath([i,j],[Modl(i-1,x),j],m)
            end
            if (posInit==l[end])
                return l
            else
                return []
            end
        elseif (m[i,j]==m[i,Modl(j+1,y)] && m[i,j]==m[i,Modl(j-1,y)])
            m[i,j]=2
            r=rand()
            if r<0.5
                l=WalkPath([i,j],[i,Modl(j+1,y)],m)
            else
                l=WalkPath([i,j],[i,Modl(j-1,y)],m)
            end
            if (posInit==l[end])
                return l
            else
                return []
            end
        elseif (m[i,j]==m[Modl(i+1,x),j] && m[i,j]==m[i,Modl(j+1,y)])
            m[i,j]=2
            r=rand()
            if r<0.5
                l=WalkPath([i,j],[Modl(i+1,x),j],m)
            else
                l=WalkPath([i,j],[i,Modl(j+1,y)],m)
            end
            if (posInit==l[end])
                return l
            else
                return []
            end
        elseif (m[i,j]==m[Modl(i+1,x),j] && m[i,j]==m[i,Modl(j-1,y)])
            m[i,j]=2
            r=rand()
            if r<0.5
                l=WalkPath([i,j],[Modl(i+1,x),j],m)
            else
                l=WalkPath([i,j],[i,Modl(j-1,y)],m)
            end
            if (posInit==l[end])
                return l
            else
                return []
            end
        elseif (m[i,j]==m[Modl(i-1,x),j] && m[i,j]==m[i,Modl(j+1,y)])
            m[i,j]=2
            r=rand()
            if r<0.5
                l=WalkPath([i,j],[Modl(i-1,x),j],m)
            else
                l=WalkPath([i,j],[i,Modl(j+1,y)],m)
            end
            if (posInit==l[end])
                return l
            else
                return []
            end
        elseif (m[i,j]==m[Modl(i-1,x),j] && m[i,j]==m[i,Modl(j-1,y)])
            m[i,j]=2
            r=rand()
            if r<0.5
                l=WalkPath([i,j],[Modl(i-1,x),j],m)
            else
                l=WalkPath([i,j],[i,Modl(j-1,y)],m)
            end
            if (posInit==l[end])
                return l
            else
                return []
            end
        else
            return []
        end
    end
    function SearchAllCycles1(m;test=false)
        x=size(m)[1]
        y=size(m)[2]
        vis=zeros(x,y)
        m=RemoveTails(m)
        fin=[]
        while (sum(vis)<x*y)
            i=rand(1:x)
            j=rand(1:y)
            if test
                print("rand pos: ")
                println([i,j])
                println("vis")
                println(vis)
            end
            if vis[i,j]==0
                l=[]
                try l=CheckCycle([i,j],m)
                catch StackOverflowError
                    l=[]
                end
                if isempty(l)
                    vis[i,j]=1
                else
                    if length(l)>2
                        append!(fin,l)
                    end
                    for pos in l 
                        if test
                            print("pos: ")
                            println(pos)
                            print("tail: ")
                            println(l)
                            println("vis")
                            println(vis)
                            println("clean")
                            println(fin)
                        end
                        vis[pos[1],pos[2]]=1
                    end
                end
            end
        end
        return Set(fin)
    end
    function WalkSimplePath(latt,pos,neigLatt;test=false) 
        if (latt[pos]==0 || latt[pos]==-1 || latt[pos]==2)
            return []
        end
        l=[]
        initPos=pos
        nm=deepcopy(neigLatt)
        neig=nm[pos]
        sameNeigs=[x for x in neig if latt[x]==latt[pos]]
        if ~isempty(sameNeigs)
            push!(l,pos)
            newPos=rand(sameNeigs)
            index=Auxiliar.Index(nm[newPos],pos)
            deleteat!(nm[newPos],index)
            while initPos!=newPos && ~in(newPos,l)
                pos=newPos
                neig=nm[pos]
                sameNeigs=[x for x in neig if latt[x]==latt[pos]]
                if ~isempty(sameNeigs)
                    push!(l,pos)
                    newPos=rand(sameNeigs)
                    index=Auxiliar.Index(nm[newPos],pos)
                    deleteat!(nm[newPos],index)
                end
            end
            return l
        else
            return []
        end
    end
    function WalkComplicatedPath(latt,pos,neigLatt;test=false) 
        #if in(pos,neigLatt[pos])
        #    nm=copy(neigLatt)
        #    for b in nm[pos]
        #        if b!=pos
        #            i=Auxiliar.Index(nm[b],pos)
        #            deleteat!(nm[b],i)
        #        end
        #    end
        #    nm[pos]=[]
        #    return ([pos],nm)
        #end
        l=WalkSimplePath(latt,pos,neigLatt;test=false)
        nm=deepcopy(neigLatt)
        if length(l)==1
            nm[pos]=[]
            return ([l[1],l[1]],nm)
        end
        if ~isempty(l)
            b=[]
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
            #=
            for p in l
                append!(b,[neig for neig in neigLatt[p] if (latt[p]==latt[neig] && ~in(neig,l))])
            end
            =#
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
    end
    function SearchCycles(latt,neigLatt)
        s=size(latt)
        vis=zeros(s)
        nm=deepcopy(neigLatt)
        fin=[]
        x=[]
        while sum(vis) < prod(s)
            pos=Auxiliar.RandomPosition(latt)
            if vis[pos]!=1
        
        
                x=WalkComplicatedPath(latt,pos,nm)
        

                l=x[1]
                nm=deepcopy(x[2])
                for p in l
                    vis[p]=1
                end
                if length(l)>=2
                    append!(fin,l)
                end
            end
        end
        return (Set(fin),x[2])
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
    =#
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
    m3=RemoveTails(m3)
    neigLatt=Auxiliar.NeighborIndexLattice(m3,Auxiliar.SquareLatticeNeighborsIndex)
    pos=CartesianIndex((3,5))

    @time x=SearchCycles(m3,neigLatt)
    println(x[1])
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