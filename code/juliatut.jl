using Printf
using Statistics

todolist=3

function reminder(x)
    if x>0
        println("deje de huevonear")
    elseif x==0
        println("tienes chance de huevonear")
    end
end

reminder(todolist)
