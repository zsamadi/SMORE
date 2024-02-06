function [node_list, length_list]=radial_walk(gStruct, N,P)
% random walk where each point is progressively farther from the starting point

%     """This function returns sampled paths from a graph.
%     It makes sure points get increasing farther from
%     the starting point.
%     Input:
%     G: the graph
%     N: number of iterations
%     P: probability of ending the walk at each step (%)"""
        XYCoords=gStruct.coordinates;

    node_list = cell(N, 1);
    word_list = cell(N, 1);
    length_list= zeros(N, 1);
    m=1;

    while m<=N

                %pick a random starting node
        start = randi(length(gStruct.neighs));

        %add the index and cell type of the starting node to lists
        temp_node=zeros(1,ceil(1/P^2));
        temp_node(1) = start;
        temp_dist = 0;

        %check if the node has any neighbors
        neighbors = gStruct.neighs{start};
        iWlk=2;


        if ~isempty(neighbors) 

            while rand>P 
                neighbors =  gStruct.neighs{start};
                w = gStruct.weights{start};
                probs=w/sum(w);
                next_nodei =randsample(length(probs),1,true,probs) ;
                next_node=neighbors(next_nodei);
                next_dist = norm(XYCoords(next_node, :)-XYCoords(temp_node(1), :));

                if next_dist > temp_dist
                    temp_node(iWlk)=next_node;
                    start = next_node;
                    temp_dist = next_dist;
                    iWlk=iWlk+1;
                end
            end
        end

% a zero is added as a separator
       temp_node=temp_node(1:iWlk-1);
       if iWlk>7
            node_list{m}=temp_node;
             word_list{m}=cTypes(temp_node(1:iWlk-1));
             length_list(m)=length(temp_node);
             m=m+1;
       end
      
    end


