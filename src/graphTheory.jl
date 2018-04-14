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
    m5=rand([-1,1],100,100)
    @time x=SearchAllCycles1(m5,test=false)
    println(x)
    =#
    #=
    g=Erdos.Graph(ConstructAdjacencyMatrix(m))
    x=Erdos.conected_components(g)
    =#