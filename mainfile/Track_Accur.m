%calculates spine accuracy by using history matrix coming from
%association results.

function [track_accur,lost_tp,first_appear,allAppearMtx] = Track_Accur(history_all,tp,ns)

[m,n] = size(history_all);
track_accur = zeros(m,1);
allAppearMtx = zeros(m,n);
dummy =[];
dummy2 = [];
for r = 1:m
    count = 0;
    for c = 1:tp
        if history_all(r,c) == 1
            count = count+1;
            dummy2 = [dummy2,c];
            allAppearMtx(r,c) = c;
        else if history_all(r,c) == -1
                dummy = [dummy,c];
                allAppearMtx(r,c) = 0;
            end    
        end
    end
    track_accur(r,:) = count;
    lost_tp(r).h = dummy;
    first_appear(r).h = dummy2(1);
    dummy = [];
    dummy2 = [];
end

track_accur = track_accur/tp;

