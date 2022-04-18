function trackSpeed(r,r_start,r_end,dt,dx,h,counter)
speed = 0
counter_1 = 1
counter_2 = 1

    for j = 1:counter-1
        if r[r_start,j] < h
            counter_1 += 1
        else
            break

        end
    end

    counter_2 = counter_1
    for i = counter_1:counter-1

        if r[r_end,i] < h

        counter_2 += 1

        else
            break


        end


    speed = (((r_end-r_start))/((counter_2 - counter_1)))*(dr/dt)
end

    return speed

end
