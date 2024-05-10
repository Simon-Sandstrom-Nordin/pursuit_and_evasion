function fieldmap = getFieldmap(spelplan, itterationer)
%tar som indata spelplanen med endast väggar utsatta. Väggar antas representeras av 1:or 
% och resterande pukter som 0:or. Tar även antalet itterationer vi önskar avända. 
storlek = size(spelplan, 1); %tar storleken på spelplanen 

for i = 1:itterationer %antal itteratner vi gör
    for x = 2:storlek-1 %antas vara väggar runtomkring planen 
        for y = 2:storlek-1 %antas vara väggar runtomkring planen 
            if spelplan(x,y) == 1 %För att skippa väggarna. Dvs hålla deras pot konstant
                continue
            end

            %här ändrar vi värdet, konstanten vi delar med är oviktig.
            %Krävs ändast att den är större än 4. Annars konvergerar det
            %mot samma potential överallt eller mot oändligheten. Så för
            %gradientens skull, vilket är det ända vi bryr oss om, så är
            %spelar det ingen roll. Blir bara en parameter för hur mycket
            %vi bryr oss om att hålla oss borta från väggar/hur bekväma vi
            %är med att komma nära väggar
            spelplan(x,y) = (spelplan(x-1,y)+spelplan(x+1,y)+spelplan(x,y-1)+spelplan(x,y+1))/4.05; 

        end
    end
end


% fieldmap_log = spelplan;
% 
% for k = 1:length(fieldmap_log)
%     for j = 1:length(fieldmap_log)
%         fieldmap_log(k,j) = fieldmap_log(k,j);
%     end
% end
% 
% writematrix(fieldmap_log,'potential_map.csv')

fieldmap = spelplan;
writematrix(fieldmap,"map.csv")
end