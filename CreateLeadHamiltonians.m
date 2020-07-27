%%   Eotvos Quantum Transport Utilities - CreateLeadHamiltonians
%    Copyright (C) 2016 Peter Rakyta, Ph.D.
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see http://www.gnu.org/licenses/.
%
%> @addtogroup basic Basic Functionalities
%> @{
%> @file CreateLeadHamiltonians.m
%> @brief Class to create and store Hamiltonian of the translational invariant leads.
%> @} 
%> @brief Class to create and store Hamiltonian of the translational invariant leads.
%> The notations and the structure of the Hamiltonian is defined accroding to the following image:
%> @image html Lead_Hamiltonian.jpg
%> @image latex Lead_Hamiltonian.jpg
%> @Available
%> EQuUs v4.8 or later
%%
classdef CreateLeadHamiltonians < handle & Messages
    
    
    properties ( Access = protected )
    %> An instance of the structure param
    param
    %> The orientation of the lead. Set +1 is the "incoming" direction of the propagating states is defined in the +x or +y direction, and "-1" otherwise.
    Lead_Orientation
    %> The id number of the current lead.
    Hanyadik_Lead
    %> The number of the sites in the cross section.
    M 
    %> An instance of the structure lead_param
    params
    %> The tranverse momentum for transverse computations.
    q
    %> List of sites in the unit cell that should be kept after decimation.
    kulso_szabfokok
    %> An instance of the structure coordinates.
    coordinates
    %> K0=H0-E*S0, see Eq (4) of PRB 78, 035407 
    K0
    %> K1=H1-E*S1, see Eq (4) of PRB 78, 035407 
    K1
    %> K1adj=H1adj-E*S1', see Eq (4) of PRB 78, 035407 
    K1adj
    %> K1_transverse=H1_transverse-E*S1_transverse
    K1_transverse   
    %> K1_skew_right = H1_skew_right - E*S1_skew_right.
    K1_skew_right
    %> K1_skew_left = H1_skew_left - E*S1_skew_left.
    K1_skew_left    
    %> The Hamiltonian of a unit cell.
    H0
    %> The coupling Hamiltonian between the unit cells.
    H1
    %> The coupling Hamiltonian between the unit cells in the opposite direction as H1. (For complex energies they differ from each other.)
    H1adj
    %> Obsolete
    H00 
    %> The transverse coupling between the slabs for transverse calculations.
    H1_transverse
    %> The skew coupling in the right direction between Hamiltonians H0 for transverse calculations.
    H1_skew_right
    %> The skew coupling in the left direction between Hamiltonians H0 for transverse calculations.
    H1_skew_left
    %> The overlap integrals of a unit cell.
    S0
    %> The overlap integrals between the unit cells.
    S1
    %> The adjungate of the overlap integrals between the unit cells.
    S1adj    
    %> The overlap integrals between the slabs for transverse calculations.
    S1_transverse
    %> The overlap integrals for the skew coupling in the right direction between Hamiltonians H0 for transverse calculations.
    S1_skew_right
    %> The overlap integrals for the skew coupling in the left direction between Hamiltonians H0 for transverse calculations.
    S1_skew_left    
    %> The matrix of the Peierls phases in the unit cell.
    fazis_mtx_H0
    %> The matrix of the Peierls phases in the coupling matrix between the unit cells.
    fazis_mtx_H1
    %> The matrix of the Peierls phases in the transverse coupling matrix between the unit cells.
    fazis_mtx_H1t
    %> The matrix of the Peierls phases in the skew coupling matrix in the right direction between the unit cells.
    fazis_mtx_H1_skew_right
    %> The matrix of the Peierls phases in the skew coupling matrix in the left direction between the unit cells.
    fazis_mtx_H1_skew_left    
    %> A logical value. True if the Hamiltonians were created, false otherwise.
    HamiltoniansCreated
    %> A logical value. True if the Hamiltonians were decimated, false otherwise.
    HamiltoniansDecimated 
    %> A logical value. True if the overlap integrals were applied, false otherwise.
    OverlapApplied    
    %> A logical value. True if magnetic field was applied in the Hamiltonians, false otherwise.
    MagneticFieldApplied
    %> A logical value. True if a gauge transformation was incorporated into the Hamiltonians or false otherwise.
    GaugeTransformationApplied
    %> list of optional parameters (see http://www.mathworks.com/help/matlab/ref/varargin.html for details)
    varargin
    end


methods ( Access = public )
%% constructorof the class
%> @brief Constructor of the class.
%> @param Opt An instance of the structure Opt.
%> @param param An instance of structure param.
%> @param varargin Cell array of optional parameters. See #InputParsing for details.
%> @return An instance of the class
    function obj = CreateLeadHamiltonians(Opt, param, varargin)    
        obj = obj@Messages( Opt );
        obj.param = param;
        obj.varargin = varargin;

        
	    obj.Initialize();
        
    end
    
    
    
%% ApplyOverlapMatrices
%> @brief Applies the overlap matrices to the Hamiltonians: K = H-ES
%> @param E The energy value.
    function ApplyOverlapMatrices(obj, E)
        
        if obj.OverlapApplied
            obj.display('EQuUs:CreateLeadHamiltonians:ApplyOverlapMatrices: Overlap matrices were already applied.');
            return;
        end
        
        if ~isempty( obj.S0 )
            obj.K0 = obj.H0 - E*obj.S0;     
            
        else
            obj.K0 = obj.H0 - speye(size(obj.H0))*E;
        end
        
        if ~isempty( obj.S1 )
            obj.K1 = obj.H1 - E*obj.S1;
            obj.K1adj = obj.H1adj - E*obj.S1';
        else
            obj.K1 = obj.H1;
            obj.K1adj = obj.H1adj;
        end
        
        
        if ~isempty( obj.S1_transverse )
            obj.K1_transverse = obj.H1_transverse - E*obj.S1_transverse;
        elseif ~isempty( obj.H1_transverse )
            obj.K1_transverse = obj.H1_transverse;
        end  
        
        
        if ~isempty( obj.S1_skew_right )
            obj.K1_skew_right = obj.H1_skew_right - E*obj.S1_skew_right;
        elseif ~isempty( obj.H1_skew_right )
            obj.K1_skew_right = obj.H1_skew_right;
        end      
        
        
        if ~isempty( obj.S1_skew_left )
            obj.K1_skew_left = obj.H1_skew_left - E*obj.S1_skew_left;
        elseif ~isempty( obj.H1_skew_left )
            obj.K1_skew_left = obj.H1_skew_left;
        end      
        
        
        obj.OverlapApplied = true;                
    end
    
    
%% CreateHamiltonians
%> @brief Creates the Hamiltonians H_0 and H_1 of the lead. The created Hamiltonians are stored by within the object. 
%> @param varargin Cell array of optional parameters (https://www.mathworks.com/help/matlab/ref/varargin.html):
%> @param 'toSave' Logical value. If true, the created Hamiltonians are saved into a file 'Hamiltoni_Lead_' + num2str(Hanyadik_Lead) + '.mat'.
%> @param 'CustomHamiltonian' An instance of class #Custom_Hamiltonians describing external source of Hamiltonians. 
    function CreateHamiltonians(obj, varargin )
        
        p = inputParser;
        p.addParameter('toSave', 0);
        p.addParameter('CustomHamiltonian', []);  
        p.parse(varargin{:});

        toSave  = p.Results.toSave;
        CustomHamiltonian = p.Results.CustomHamiltonian;
        
        obj.setM();
        
        if ~isempty( obj.Opt.custom_Hamiltonians )
            if isempty( CustomHamiltonian )
                CustomHamiltonian = Custom_Hamiltonians( obj.Opt, obj.param ); 
            end
        
            if ~CustomHamiltonian.Read( 'Hamiltonians_loaded' )
                CustomHamiltonian.LoadHamiltonians();
            end
            
            obj.coordinates   = CustomHamiltonian.Read( 'coordinates' );
            obj.coordinates   = obj.coordinates{obj.Hanyadik_Lead};
            obj.H0            = CustomHamiltonian.Read( 'H0' );
            obj.H0            = obj.H0{obj.Hanyadik_Lead};
            obj.H1            = CustomHamiltonian.Read( 'H1' );
            obj.H1            = obj.H1{obj.Hanyadik_Lead};
            obj.H1adj         = obj.H1';
            obj.H1_transverse = CustomHamiltonian.Read( 'H1_transverse' ); 
            obj.H1_transverse = obj.H1_transverse{obj.Hanyadik_Lead};
            obj.S0            = CustomHamiltonian.Read( 'S0' );
            if iscell( obj.S0 )
                obj.S0            = obj.S0{obj.Hanyadik_Lead};
            end
            obj.S1            = CustomHamiltonian.Read( 'S1' );
            if iscell( obj.S1 )
                obj.S1            = obj.S1{obj.Hanyadik_Lead};
                obj.S1adj         = obj.S1';
            end
            obj.S1_transverse = CustomHamiltonian.Read( 'S1_transverse' ); 
            if iscell( obj.S1_transverse )
                obj.S1_transverse = obj.S1_transverse{obj.Hanyadik_Lead};      
            end
            obj.M             = size(obj.H0,1);
            %obj.kulso_szabfokok = 1:obj.M;
        elseif strcmpi(obj.Opt.Lattice_Type, 'Square')
            createH = Square_Lead_Hamiltonians();
            [obj.H0,obj.H1,obj.H1_transverse,obj.coordinates] = createH.Square_Hamiltonians(obj.params, obj.M );
            obj.H1adj = obj.H1';
            obj.kulso_szabfokok = 1:obj.M;
        elseif strcmpi(obj.Opt.Lattice_Type, 'SSH')
            createH = Square_Lead_Hamiltonians();
            [obj.H0,obj.H1,obj.H1_transverse,obj.coordinates] = createH.SSH_Hamiltonians(obj.params );
            obj.H1adj = obj.H1';
            obj.kulso_szabfokok = 1;
        elseif strcmp(obj.Opt.Lattice_Type, 'Lieb')
            createH = Square_Lead_Hamiltonians();
            [obj.H0,obj.H1,obj.H1_transverse,obj.coordinates] = createH.Lieb_Hamiltonians(obj.params, obj.M );
            obj.H1adj = obj.H1';
            obj.M = size(obj.H0,1);
            obj.kulso_szabfokok = [];%1:3:3*obj.M-2;    
        elseif strcmp(obj.Opt.Lattice_Type, 'BiTeI')
            createH = BiTeI_Lead_Hamiltonians();
            [obj.H0,obj.H1,obj.H1_transverse,obj.coordinates] = createH.BiTeI_Hamiltonians(obj.params, obj.M );
            obj.H1adj = obj.H1';
            obj.M = size(obj.H0,1);
            obj.kulso_szabfokok = [];%sort([1:4:2*obj.M-1, 2:4:2*obj.M], 'ascend');%[];
        elseif strcmp(obj.Opt.Lattice_Type, 'Graphene') || strcmpi(obj.Opt.Lattice_Type, 'H')
            createH = Hex_Lead_Hamiltonians();
            [obj.H0,obj.H1,obj.H1_transverse,obj.coordinates] = createH.Graphene_Hamiltonians(obj.params, obj.M, obj.params.End_Type );
            obj.H1adj = obj.H1';
            obj.kulso_szabfokok = (1:obj.M);
            
            %{a
            x=obj.coordinates.x;

            index_eliminate =  sort(find(x == max(x)),'descend');
            
            obj.coordinates.x(index_eliminate) = [];
            obj.coordinates.y(index_eliminate) = [];
            
            for i=1:length(index_eliminate)
               obj.H0(index_eliminate(i),:)=[]; 
               obj.H0(:,index_eliminate(i))=[];
               
               obj.H1(index_eliminate(i),:)=[]; 
               obj.H1(:,index_eliminate(i))=[];
               
               obj.H1_transverse(index_eliminate(i),:)=[]; 
               obj.H1_transverse(:,index_eliminate(i))=[];
            end  
            obj.kulso_szabfokok = (1:obj.M-1);
            %}
            %{
            obj.H0 = H0_tmp;
            obj.H1 = H1_tmp;
            obj.H1_transverse = H1_transverse_tmp;
            obj.coordinates.x = 
            
                        
            value = -1e-5;

            [ii,jj,v] = find(obj.H0);

            for k = 1:length(ii)
                for kk= 1:length(index_eliminate)
                    if( ii(k) == index_eliminate(kk) || jj(k) == index_eliminate(kk))
                        obj.H0(ii(k),jj(k))=value;
                    end
                end
            end

            [ii,jj,v] = find(obj.H1);

            for k = 1:length(ii)
                for kk= 1:length(index_eliminate)
                    if( ii(k) == index_eliminate(kk) || jj(k) == index_eliminate(kk))
                        obj.H1(ii(k),jj(k))=value;
                    end
                end
            end

            [ii,jj,v] = find(obj.H1_transverse);

            for k = 1:length(ii)
                for kk= 1:length(index_eliminate)
                    if( ii(k) == index_eliminate(kk) || jj(k) == index_eliminate(kk))
                        obj.H1_transverse(ii(k),jj(k))=value;
                    end
                end
            end
            %}    
        elseif strcmp(obj.Opt.Lattice_Type, 'Graphene_SOC')
            createH = Graphene_SOC_Lead_Hamiltonians();
            [obj.H0, obj.H1, obj.H1_transverse, obj.H1_skew_left, obj.H1_skew_right, obj.coordinates] = createH.Graphene_SOC_Hamiltonians(obj.params, obj.M, obj.params.End_Type );
            obj.H1adj = obj.H1';
            obj.kulso_szabfokok = [1:2*obj.M, size(obj.H0,1)/2 + (1:2*obj.M)];
        elseif strcmpi(obj.Opt.Lattice_Type, 'Graphene_Bilayer')
            createH = Graphene_Bilayer_Lead_Hamiltonians();
            [obj.H0,obj.H1,obj.H1_transverse,obj.coordinates] = createH.Graphene_Bilayer_Hamiltonians(obj.params, obj.M, obj.params.End_Type );
            obj.H1adj = obj.H1';
            obj.kulso_szabfokok = [1:obj.M, size(obj.H0,1)/2 + (1:obj.M)];
            obj.M = 2*obj.M;
        elseif strcmpi(obj.Opt.Lattice_Type, 'Graphene_Bilayer_2')
            createH = Graphene_Bilayer_Lead_Hamiltonians();
            [obj.H0,obj.H1,obj.H1_transverse,obj.coordinates] = createH.Graphene_Bilayer_Hamiltonians2(obj.params, obj.M, obj.params.End_Type );
            obj.H1adj = obj.H1';
            obj.kulso_szabfokok = [1:obj.M, size(obj.H0,1)/2 + (1:obj.M)];
            obj.M = 2*obj.M;         
        elseif strcmpi(obj.Opt.Lattice_Type, 'Graphene_Bilayer_3')
            createH = Graphene_Bilayer_Lead_Hamiltonians();
            [obj.H0,obj.H1,obj.H1_transverse,obj.coordinates] = createH.Graphene_Bilayer_Hamiltonians3(obj.params, obj.M, obj.params.End_Type );
            obj.H1adj = obj.H1';
            obj.kulso_szabfokok = [1:obj.M+1, size(obj.H0,1)/2 + (1:obj.M+1)];
            obj.M = 2*(obj.M+1);           
        elseif strcmpi(obj.Opt.Lattice_Type, 'Silicene')
            createH = Silicene_Lead_Hamiltonians();
            [obj.H0,obj.H1,obj.H1_transverse,obj.coordinates] = createH.Silicene_Hamiltonians(obj.params, obj.M, obj.params.End_Type );
            obj.H1adj = obj.H1';
            obj.kulso_szabfokok = [1:obj.M, size(obj.H0,1)/2 + (1:obj.M)];          
            obj.M = 2*obj.M;
        elseif strcmpi(obj.Opt.Lattice_Type, 'Triangle')
            createH = Triangle_Lead_Hamiltonians();
            [obj.H0, obj.H1 ,obj.H1_transverse, obj.H1_skew_left, obj.coordinates] = createH.Triangle_Hamiltonians(obj.params, obj.M );
            obj.H1adj = obj.H1';
            obj.kulso_szabfokok = [1:obj.M];  
        elseif strcmpi(obj.Opt.Lattice_Type, 'TMDC_Monolayer')
            createH = TMDC_Monolayer_Lead_Hamiltonians();
            [obj.H0, obj.H1 ,obj.H1_transverse, obj.H1_skew_left, obj.coordinates] = createH.TMDC_Monolayer_Hamiltonians(obj.params, obj.M );
            obj.H1adj = obj.H1';
            obj.kulso_szabfokok = 1:size(obj.H0,1);
            obj.M = size(obj.H0,1);
        elseif strcmpi(obj.Opt.Lattice_Type, 'TMDC_Monolayer_SOC')
            createH = TMDC_Monolayer_SOC_Lead_Hamiltonians();
            [obj.H0, obj.H1 ,obj.H1_transverse, obj.H1_skew_left, obj.coordinates] = createH.TMDC_Monolayer_SOC_Hamiltonians(obj.params, obj.M );
            obj.H1adj = obj.H1';
            obj.kulso_szabfokok = 1:size(obj.H0,1);
            obj.M = size(obj.H0,1);
        elseif strcmpi(obj.Opt.Lattice_Type, 'TMDC_Bilayer_SOC')
            createH = TMDC_Bilayer_SOC_Lead_Hamiltonians();
            [obj.H0, obj.H1 ,obj.H1_transverse, obj.H1_skew_left, obj.coordinates] = createH.TMDC_Bilayer_SOC_Hamiltonians(obj.params, obj.M );
            obj.H1adj = obj.H1';
            obj.kulso_szabfokok = 1:size(obj.H0,1);
            obj.M = size(obj.H0,1);
        else
            error(['EQuUs:', class(obj), ':CreateHamiltonians'], 'Unrecognized lattice type, or a valid custom source for the Hamiltonians was not set.')
        end         
            
        obj.Transform2Spin();
        obj.Transform2BdG();
                        
        if toSave
            saveLeads();
        end
            
        obj.HamiltoniansCreated = true;
        obj.OverlapApplied = false;
        obj.HamiltoniansDecimated = false;
        obj.MagneticFieldApplied  = false;
        obj.GaugeTransformationApplied = false; 
        
        
    end
    
%% Transform2Spin
%> @brief Transforms the Hamiltonians and the overlap matrices to include electron spin.
    function Transform2Spin( obj )
    
        if isempty(obj.Opt.Spin) || ~obj.Opt.Spin
            return
        end
        
        if ~isempty( obj.coordinates.spinup )
            % spin is also included in the model
            return
        end
        
        fnames = fieldnames( obj.coordinates );
        for idx = 1:length(fnames)
            fname = fnames{idx};
            if strcmp(fname, 'a') || strcmp(fname, 'b')
                continue
            end                    
            obj.coordinates.(fname) = [ obj.coordinates.(fname); obj.coordinates.(fname) ];
        end
                
        obj.coordinates.spinup = [ true(size(obj.H0,1),1); false(size(obj.H0,1),1)];
                
        obj.kulso_szabfokok = [obj.kulso_szabfokok, obj.kulso_szabfokok+size(obj.H0,1)];
               
        
        % transforming the Hamiltonians        
        obj.H0 = [obj.H0, sparse([],[],[], size(obj.H0,1), size(obj.H0,2)); sparse([],[],[], size(obj.H0,1), size(obj.H0,2)), obj.H0];
        obj.H1 = [obj.H1, sparse([],[],[], size(obj.H1,1), size(obj.H1,2)); sparse([],[],[], size(obj.H1,1), size(obj.H1,2)), obj.H1];
        obj.H1adj = [obj.H1adj, sparse([],[],[], size(obj.H1adj,1), size(obj.H1adj,2)); sparse([],[],[], size(obj.H1adj,1), size(obj.H1adj,2)), obj.H1adj];
        
        obj.H1_transverse = [obj.H1_transverse, sparse([],[],[], size(obj.H1_transverse,1), size(obj.H1_transverse,2)); ...
                sparse([],[],[], size(obj.H1_transverse,1), size(obj.H1_transverse,2)), obj.H1_transverse];
        
        % transforming th eoverlap integrals  
        if ~isempty( obj.S0 )
            obj.S0 = [obj.S0, sparse([], [], [], size(obj.S0,1), size(obj.S0,2)); ...
                sparse([], [], [], size(obj.S0,1), size(obj.S0,2)), obj.S0];
        end
        
        if ~isempty(obj.S1)
            obj.S1 = [obj.S1, sparse([], [], [], size(obj.S1,1), size(obj.S1,2)); ...
                sparse([], [], [], size(obj.S1,1), size(obj.S1,2)), obj.S1];
        end
        
        if ~isempty( obj.S1_transverse )
            obj.S1_transverse = [obj.S1_transverse, sparse([], [], [], size(obj.S1_transverse,1), size(obj.S1_transverse,2)); ...
                sparse([], [], [], size(obj.S1_transverse,2), size(obj.S1_transverse,1)), obj.S1_transverse];
        end
                
                
        obj.M = 2*obj.M;        
        
    end    
    
%% Transform2BdG
%> @brief Transforms the Hamiltonians and the overlap matrices into the BdG model in the Nambu space representation according to 
%> <a href="http://iopscience.iop.org/article/10.1088/1367-2630/9/8/278/meta">New Journal of Physics 9 (2007) 278</a>. 
%> It is assumed, that the Hamiltonian is already transfromed to the grand canonical operator: \f$ \hat{H} \rightarrow \hat{H} - E_F\hat{N}\f$
    function Transform2BdG( obj )
    
        if isempty(obj.Opt.BdG) || ~obj.Opt.BdG
            obj.display(['EQuUs:', class(obj), ':Transform2BdG: BdG option is not set to true in the computational parameters.'])
            return
        end
        
        if ~isempty( obj.coordinates.BdG_u )
            % already transformed into BdG model
            obj.display(['EQuUs:', class(obj), ':Transform2BdG: Hamiltonians already transformed intt the BdG model.'])
            return
        end
                
        % transforming the coordinates
        obj.coordinates = obj.coordinates.Transform2BdG();
        % checking the crated BdG array ---- this is necessary when the structure coordinates was nat filled in with informations, when the Hamiltonians were created
        if isempty( obj.coordinates.BdG_u )
            obj.coordinates.BdG_u = [ true(size(obj.H0,1),1); false(size(obj.H0,1),1)];
        end
                
        obj.kulso_szabfokok = [obj.kulso_szabfokok, obj.kulso_szabfokok+size(obj.H0,1)];
               
        pair_potential = obj.params.pair_potential;
        
        % transforming the Hamiltonians
        if isempty(obj.S0)
            S0 = speye(size(obj.H0));
        else
            S0 = obj.S0;
        end
        
        obj.H0 = [obj.H0, S0*pair_potential; S0*conj(pair_potential), -conj(obj.H0)];
        
        if isempty(obj.S1)
            S1 = sparse([],[],[], size(obj.H1,1), size(obj.H1,2));
        else
            S1 = obj.S1;
        end        
        
        obj.H1 = [obj.H1, S1*pair_potential; S1*conj(pair_potential), -conj(obj.H1)];
        obj.H1adj = [obj.H1adj, S1'*pair_potential; S1'*conj(pair_potential), -conj(obj.H1adj)];
        
        if isempty(obj.S1_transverse)
            S1_transverse = sparse([],[],[], size(obj.H1_transverse,1), size(obj.H1_transverse,2));
        else
            S1_transverse = obj.S1_transverse;
        end                
        
        obj.H1_transverse = [obj.H1_transverse, S1_transverse*pair_potential; S1_transverse*conj(pair_potential), -conj(obj.H1_transverse)];
        
        % transforming the overlap integrals  
        if ~isempty( obj.S0 )
            obj.S0 = [obj.S0, sparse([], [], [], size(obj.S0,1), size(obj.S0,2)); sparse([], [], [], size(obj.S0,1), size(obj.S0,2)), obj.S0];
        end
        
        if ~isempty(obj.S1)
            obj.S1     = [obj.S1, sparse([], [], [], size(obj.S1,1), size(obj.S1,2)); sparse([], [], [], size(obj.S1,1), size(obj.S1,2)), obj.S1];
            obj.S1adj  = obj.S1';
        end
        
        if ~isempty( obj.S1_transverse )
            obj.S1_transverse = [obj.S1_transverse, sparse([], [], [], size(obj.S1_transverse,1), size(obj.S1_transverse,2)); ...
                sparse([], [], [], size(obj.S1_transverse,2), size(obj.S1_transverse,1)), obj.S1_transverse];
        end
                
        % the number of sites in the cross section becomes twice as many as in the normal case
        obj.M = 2*obj.M;        
        
    end
        
%% CalcSpektrum
%> @brief Calculates the band structure of the lead.
%> @param varargin Cell array of optional parameters:
%> @param 'toPlot' Set 1 in order to plot the calculated spectrum, 0 (default) otherwise
%> @param 'ka_min' The lower bound of the wave numbers. (Default is -pi.)
%> @param 'ka_max' The upper bound of the wave numbers. (Default is pi.)
%> @param 'ka_num' The number of wave number points involved in the calculations. (Default is 300.)
%> @param 'ka_vec' One dimensional array of the k-pints. (Overrides previous attributes related to the k-vector array.)
%> @param 'center' The calculated energy eigenvalues are centered around this value. (Default is 0.001.)
%> @param 'db' The number of the calculated eigenvalues.
%> @param 'offset' Offset value to shift the spectrum along the energy axis.
%> @param 'calcWaveFnc' Logical value. Set true to calculate also the wave functions, or false (default) otherwise.
%> @return [1] ka_num x 2 array of the calculated spactrum. In the first column are the k-points, whil ein the second columns are the calculated energy points.
%> @return [2] The calculated wave functions stored in a structure #WaveFnc.
    function [spectrum, WaveFnc] = CalcSpektrum( obj, varargin )
        
        p = inputParser;
        p.addParameter('toPlot', 0, @isscalar);  
        p.addParameter('ka_min', -pi, @isscalar);  
        p.addParameter('ka_max', pi, @isscalar);  
        p.addParameter('ka_num', 300, @isscalar);
        p.addParameter('ka_vec', [] );
        p.addParameter('center', 0.001);  
        p.addParameter('db', min([10,size(obj.H0,1)]), @isscalar); 
		p.addParameter('offset', mean(diag(obj.H0)) ); %Offset value to shift the spectrum along the energy axis
        p.addParameter('calcWaveFnc', false ); %Logical value. Set true to calculate also the wave functions, or false (default) otherwise.
        
        p.parse(varargin{:});    
        toPlot     = p.Results.toPlot;
        ka_min     = p.Results.ka_min;
        ka_max     = p.Results.ka_max;
        ka_num     = p.Results.ka_num;
        ka_vec     = p.Results.ka_vec;
        center = p.Results.center;
        db             = p.Results.db;
		offset     = p.Results.offset;
        calcWaveFnc = p.Results.calcWaveFnc;
        
    
    % check wheter Hamiltonians are decimated or not
	if obj.HamiltoniansDecimated
		obj.display(['EQuUs:', class(obj), ':CalcSpektrum: Hamiltonians are decimated. Please recreate Hamiltonians to calculate spectrum.'], 1)
		spectrum = NaN;
		return
    end
    
    % check the number of eigenvalues
    db = min([ db, size(obj.H0,1)]);    
    
       
    obj.display('Calculating spectrum')   
    
    % discrete increment of the wavenumber array
    deltak = (ka_max-ka_min)/ka_num;
    
    % allocating temporary matricest;
    S0_loc = obj.S0;
    S1_loc = obj.S1;
    S1_transverse_loc = obj.S1_transverse;    
    
    tic    
    % creating the one-dimensional array for the wave numbers if not given
    if isempty(ka_vec)
        ka_vec = ka_min:deltak:ka_max;
    end    
    
    % allocating arrays for the results
    spectrum = cell(length(ka_vec),1);
    if calcWaveFnc
        WaveFnc = structures('WaveFnc');
    else
        WaveFnc = [];
    end
    
    % calculating the spectrum
    for idx=1:length(ka_vec)
        % obtaining the k-dependent effective Hamiltonian
        H    = obj.MomentumDependentHamiltonian( ka_vec(idx), obj.q );
        H    = H - eye(size(H))*offset;
        
        if isempty( S0_loc ) && isempty( S1_loc )
            if calcWaveFnc
                % calculations including the wavefunctions
                if db < size(H,1)-1
                    [WaveFnc_tmp, E]  = eigs(H, db, center);
                else
                    [WaveFnc_tmp, E]  = eig(H);
                end
                E = diag(E);
                WaveFnc(idx).WaveFnc = WaveFnc_tmp;
                WaveFnc(idx).E = E;
                WaveFnc(idx).ka = ka_vec(idx);
            else
                % only eigenvalues are calculated
                if db < size(H,1)-1
                    E    = eigs(H, db, center);
                else
                    E    = eig(full(H));
                end
            end
        else
            % calculations including the overlap matrices
            S    = secular_H( S0_loc, S1_loc, S1_transverse_loc, ka_vec(idx));
            if calcWaveFnc
                % calculations including the wavefunctions
                if db < size(H,1)-1
                    [WaveFnc_tmp, E]  = eigs(H, S, db, center);
                else
                    [WaveFnc_tmp, E]  = eig(H, S);
                end
                E = diag(E);
                WaveFnc.WaveFnc{idx} = WaveFnc_tmp;
                WaveFnc.E{idx} = E;
                WaveFnc.ka(idx) = ka_vec(idx);
            else
                % only eigenvalues are calculated
                if db < size(H,1)-1
                    E    = eigs(H, S, db, center);
                else
                    E    = eig(H, S);
                end
            end
        end
        
        spectrum{idx} = [ones(length(E),1)*ka_vec(idx), E];        
        
    end
    % reorganize the calculated data
    spectrum = cell2mat(spectrum);
    toc
 
    obj.display('Spectrum calculated')
    
    % check whether to plot the spectrum
    if ~toPlot
        return
    end
    
    % plot the calculated spectrum
    figure1 = figure();
    fontsize = 9;
    
    spectrum(:,2) = real(spectrum(:,2));
    x_lim = [min(spectrum(:,1)) max(spectrum(:,1))];
    y_lim = [min(spectrum(:,2)) max(spectrum(:,2))];

    Position = [0.5 0.58 0.33 0.4];
    axes_spectrum = axes('Parent',figure1, 'Position', Position,...
                'Visible', 'on',...
                'FontSize', fontsize,...
                'xlim', x_lim,...
                'ylim', y_lim,...                'XTick', XTick,...                'YTick', YTick,...
                'Box', 'on',...
                'FontName','Times New Roman');
        hold on; 
        
    plot(spectrum(:,1), spectrum(:,2),'.', 'MarkerSize', 4, 'Parent', axes_spectrum )  
        
    % Create xlabel
    xlabel_position = [0 0 0];
    xlabel('ka','FontSize', fontsize,'FontName','Times New Roman', 'Parent', axes_spectrum);
    xlabel_handle = get(axes_spectrum,'XLabel');  
    set(xlabel_handle, 'Position', get(xlabel_handle, 'Position') + xlabel_position);

    % Create ylabel
    ylabel_position = [0 0 0];
    ylabel('E [eV]','FontSize', fontsize,'FontName','Times New Roman', 'Parent', axes_spectrum);
    ylabel_handle = get(axes_spectrum,'YLabel'); 
    set(ylabel_handle, 'Position', get(ylabel_handle, 'Position') + ylabel_position);  

    
        % --------------------------------------
        % nested functions
        % Hamiltonian for the secular equation
        function H = secular_H(H0,H1, H1_transverse, k)
            
            if ~isempty(q) && ~obj.HamiltoniansDecimated
               H0 = H0 + H1_transverse*diag(exp(-1i*q)) + H1_transverse'*diag(exp(1i*q)); 
            end
            H = H0 + H1*exp(1i*k) + H1'*exp(-1i*k);
        
        end
        
        % end nested functions

 
    end   
    
%% saveLeads
%> @brief Save Lead Hamiltonians  into a file 'Hamiltoni_Lead_' + num2str(Hanyadik_Lead) + '.mat'.  
    function saveLeads( obj )
        save(['Hamiltoni_Lead_',num2str(obj.Hanyadik_Lead),'.mat'], 'H0', 'H1', 'kulso_szabfokok', 'Lead_Orientation', 'fazis_mtx_H0', 'fazis_mtx_H1', 'params');
    end        
    
%% ShiftCoordinates
%> @brief Shifts the coordinates of the sites by an integer multiple of the lattice vector #Coordinates.a.
%> @param shift Integer by which the coordinates are shifted.
    function ShiftCoordinates( obj, shift )    
        if isempty(shift) || shift == 0
            return
        end
        
        % shifting the coordinates along the translational invariant direction
        obj.coordinates = obj.coordinates.Shift( shift*obj.coordinates.a );
                
    end    
    
    
    
%% ShiftLead
%> @brief Shifts the on-site energies in the leads by a given energy.
%> @param Energy The enrgy value.
    function ShiftLead( obj, Energy )        
        obj.H0 = obj.H0 + sparse(1:size(obj.H0,1), 1:size(obj.H0,2),Energy,size(obj.H0,1),size(obj.H0,2)); 
        obj.params.epsilon = obj.params.epsilon + Energy;
        
        obj.K0 = [];
        obj.OverlapApplied = false;
    end

%% AddPotential
%> @brief Adds on-site potential to the Hamiltonian H0.
%> @param V The potential calculated on the sites.
    function AddPotential( obj, V )        

        if ( size(V,2) == 1 ) || ( size(V,1) == 1 )
            if length(V) == size(obj.H0,1)
                obj.H0 = obj.H0 + sparse(1:size(obj.H0,1), 1:size(obj.H0,2),V,size(obj.H0,1),size(obj.H0,2));
                obj.K0 = [];
                obj.OverlapApplied = false;
            else
                disp(' Wrong dimension of potential: 1')
                return
            end
        elseif norm(size(V) - size(obj.H0)) < 1e-6
            obj.H0 = obj.H0 + V;   
            obj.K0 = [];
            obj.OverlapApplied = false;
        else
            disp(' Wrong dimension of potential: 2')
            return
        end
        
        obj.display(' Potential added to lead')
    end

%% isSuperconducting
%> @brief Test, whether the lead is in the superconducting phase or not.
%> @return True if the lead is superconducting, false otherwise.
    function ret = isSuperconducting( obj )
        
        if ~obj.Opt.BdG
            ret = false;
            return;
        end
        
        if abs(obj.params.pair_potential) >= 1e-10
            ret = true;
        else
            ret = false;
        end
        
        return
            
    
    end
    
    
%% MomentumDependentHamiltonian
%> @brief Construct a momentum dependent (Fourier-transformed) Hamiltonian.
%> @param ka The longitudinal momentum times the lattice constant.
%> @param qb The transverse momentum times the transverse lattice constant.
%> @return Return with the momentum dependent Hamiltonian.
    function ret = MomentumDependentHamiltonian( obj, ka, qb )
        
        % check wheter Hamiltonians are decimated or not
        if obj.HamiltoniansDecimated
            obj.display(['EQuUs:', class(obj), ':MomentumDependentHamiltonian: Hamiltonians are decimated. Please recreate Hamiltonians.'], 1)
            ret = NaN;
            return
        end
        
        ret = obj.H0;
        
        ret = ret + obj.H1*diag(exp(1i*ka)) + obj.H1'*diag(exp(-1i*ka));
        
        if isempty(qb)
            return
        end
        
        
        if ~isempty(obj.H1_transverse)
            ret = ret + obj.H1_transverse*diag(exp(1i*qb)) + obj.H1_transverse'*diag(exp(-1i*qb));
        end
        
        if ~isempty(obj.H1_skew_right)
            ret = ret + obj.H1_skew_right*diag(exp(1i*(qb-ka))) + obj.H1_skew_right'*diag(exp(-1i*(qb-ka)));
        end  
        
        if ~isempty(obj.H1_skew_left)
            ret = ret + obj.H1_skew_left*diag(exp(1i*(qb+ka))) + obj.H1_skew_left'*diag(exp(-1i*(qb+ka)));
        end  
        
        
        return
            
    
    end        
    
    
%% qDependentHamiltonians
%> @brief Construct a K0 and K1 Hamiltonians (transformed to the grand canonical potential) for a given subspace of the transverse momentum number #q.
%> @return [1] the transverse momentum dependent Hamiltonian of the unit cell.
%> @return [2] the transverse momentum dependent Hamiltonian of the coupling boeween the unit cells.
%> @return [3] the transverse momentum dependent Hamiltonian of the coupling boeween the unit cells in the opposite direction.
    function [K0, K1, K1adj] = qDependentHamiltonians( obj )
        
        K0 = obj.K0;  
        K1 = obj.K1;
        
        % check wheter Hamiltonians are decimated or not
        if obj.HamiltoniansDecimated           
            K1adj = obj.K1adj;
            return
        end
        
        
        % applying transverse coupling
        if ~isempty(obj.q) && ~isempty(obj.K1_transverse)
            K0 = K0 + obj.K1_transverse*diag(exp(1i*obj.q)) + obj.K1_transverse'*diag(exp(-1i*obj.q));
        end
        
        % applying skew_right transverse coupling
        if ~isempty(obj.q) && ~isempty(obj.K1_skew_right)
            K1 = K1 + obj.K1_skew_right'*diag(exp(-1i*obj.q));
        end           
        
        % applying skew_left transverse coupling
        if ~isempty(obj.q) && ~isempty(obj.K1_skew_left)
            K1 = K1 + obj.K1_skew_left*diag(exp(1i*obj.q));
        end    
        
        K1adj = K1';
        
        
        return
            
    
    end 
    
    
    
%% qDependentOverlaps
%> @brief Construct a S0 and S1 overlap matrices for a given subspace of the transverse momentum number #q.
%> @return [1] the transverse momentum dependent overlap matrix of the unit cell.
%> @return [2] the transverse momentum dependent overlap matrix of the coupling between the unit cells.
%> @return [3] the transverse momentum dependent overlap matrix of the coupling boeween the unit cells in the opposite direction.
    function [S0, S1, S1adj] = qDependentOverlaps( obj )
        
        S0 = obj.S0;  
        S1 = obj.S1;
        
        % check wheter Hamiltonians are decimated or not
        if obj.HamiltoniansDecimated           
            S1adj = obj.S1adj;
            return
        end
        
        
        % applying transverse coupling
        if ~isempty(obj.q) && ~isempty(obj.S1_transverse)
            S0 = S0 + obj.S1_transverse*diag(exp(1i*obj.q)) + obj.S1_transverse'*diag(exp(-1i*obj.q));
        end
        
        % applying skew_right transverse coupling
        if ~isempty(obj.q) && ~isempty(obj.S1_skew_right)
            S1 = S1 + obj.S1_skew_right'*diag(exp(-1i*obj.q));
        end           
        
        % applying skew_left transverse coupling
        if ~isempty(obj.q) && ~isempty(obj.S1_skew_left)
            S1 = S1 + obj.S1_skew_left*diag(exp(1i*obj.q));
        end   
        
        S1adj = S1';
        
        
        return
            
    
    end     

%% CreateClone
%> @brief Creates a clone of the present class.
%> @return Returns with the cloned object.
%> @param varargin Cell array of optional parameters (https://www.mathworks.com/help/matlab/ref/varargin.html):
%> @param 'empty' Set true to create an empty clone, or false (default) to clone all atributes.
    function ret = CreateClone( obj, varargin )
        
        p = inputParser;
        p.addParameter('empty', false ); %Logical value. Set true for creating an empty class, or false (default) otherwise.
        
        p.parse(varargin{:});       
        empty    = p.Results.empty;
        
        ret = CreateLeadHamiltonians(obj.Opt, obj.param,...
                                            'Hanyadik_Lead', obj.Hanyadik_Lead,...
                                            'Lead_Orientation', obj.Lead_Orientation, ...
                                            'q', obj.q);               
        if empty
            return
        end                                        
                                        
        meta_data = metaclass(obj);
        
		for idx = 1:length(meta_data.PropertyList)
            prop_name = meta_data.PropertyList(idx).Name;
			ret.Write( prop_name, obj.(prop_name));  
        end      
        
    end   
    
%% Reset
%> @brief Resets all elements in the object.
    function Reset( obj )
        
        if strcmpi( class(obj), 'CreateLeadHamiltonians' )
            meta_data = metaclass(obj);
        
            for idx = 1:length(meta_data.PropertyList)
                prop_name = meta_data.PropertyList(idx).Name;
                if strcmp(prop_name, 'Opt') || strcmp( prop_name, 'param') || strcmp(prop_name, 'varargin')
                    continue
                end                
                obj.Clear( prop_name ); 
            end   
        end        
        
        obj.Initialize();
        

    end    
    

    
%% Write
%> @brief Sets the value of an attribute in the interface.
%> @param MemberName The name of the attribute to be set.
%> @param input The value to be set.
    function Write(obj, MemberName, input)
        
        if strcmp(MemberName, 'param')
            obj.param = input;
            obj.Reset()
            return
        elseif strcmpi(MemberName, 'params')
            obj.params = input;
            if isempty(obj.Hanyadik_Lead)
                obj.param.scatter = input;
            else
                obj.param.Leads{obj.Hanyadik_Lead} = input;
            end
        else
            try
        		obj.(MemberName) = input;
			catch
				error(['EQuUs:', class(obj), ':Read'], ['No property to write with name: ',  MemberName]);
			end
        end
                
    end
%% Read
%> @brief Query for the value of an attribute in the interface.
%> @param MemberName The name of the attribute to be set.
%> @return Returns with the value of the attribute.
    function ret = Read(obj, MemberName)        
        
        try
        	ret = obj.(MemberName);
		catch
			error(['EQuUs:', class(obj), ':Read'], ['No property to read with name: ',  MemberName]);
		end 
		
    end
%% Clear
%> @brief Clears the value of an attribute in the interface.
%> @param MemberName The name of the attribute to be cleared.
    function Clear(obj, MemberName)
        
        try
        	obj.(MemberName) = [];
		catch
			error(['EQuUs:', class(obj), ':Clear'], ['No property to clear with name: ',  MemberName]);
		end  
        
    end    
    
end % public methods


methods ( Access = protected )
   
%% setM
%> @brief Updates the number of sites in the cross section.
    function setM( obj )
        if isempty( obj.Hanyadik_Lead )
            obj.M             = obj.param.scatter.shape.width;
        else
            obj.M             = obj.param.Leads{obj.Hanyadik_Lead}.M;
        end  
    end
    
    
%% Initialize
%> @brief Initializes object properties.
    function Initialize(obj)
        obj.InputParsing( obj.varargin{:});

        
        obj.HamiltoniansCreated   = false; 
        obj.HamiltoniansDecimated = false;
        obj.OverlapApplied = false;
        obj.MagneticFieldApplied  = false;
        obj.GaugeTransformationApplied = false;
        
        
        
        obj.setM();
        
        if isempty( obj.Hanyadik_Lead )
            obj.params        = obj.param.scatter;  %Lead parameters
        else
            obj.params        = obj.param.Leads{obj.Hanyadik_Lead};  %Lead parameters
        end
        
    end    
    
end % protected methods
    
    
methods (Access=protected)   
    

%% InputParsing
%> @brief Parses the optional parameters for the class constructor.
%> @param varargin Cell array of optional parameters (https://www.mathworks.com/help/matlab/ref/varargin.html):
%> @param 'Hanyadik_Lead' The ID number of the current lead. Set to empty (default) for using parameters of the scatter region.
%> @param 'Lead_Orientation' Orientation of the lead. Set +1 (default) is the "incoming" direction of the propagating states is defined in the +x or +y direction, and "-1" otherwise.
%> @param 'q' The transverse momentum. Set to empty (default) for computations without transverse momentums.
    function InputParsing( obj, varargin )
        p = inputParser;
        p.addParameter('Hanyadik_Lead', []);
        p.addParameter('Lead_Orientation', 1);
        p.addParameter('q', []);
        
        % keeping unmatched attributes that possibly comes from the derived classes
        p.KeepUnmatched = true;
        
        p.parse(varargin{:});

        obj.Lead_Orientation  = p.Results.Lead_Orientation;
        obj.Hanyadik_Lead     = p.Results.Hanyadik_Lead;
        obj.q                 = p.Results.q;
        

        
        
    end

end % private methods


    
    
    
end
