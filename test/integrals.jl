using Test
using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
using Integrals
I = Integrals

@testset "Integration Tests" begin
	@testset "M" begin
		@test I.M(0,0)==0.5
		@test I.M(1,0)==0.16666666666666666
		@test I.M(2,0)==0.08333333333333333
		@test I.M(3,0)==0.05
		@test I.M(1,1)==0.041666666666666685
		@test I.M(2,0)==0.08333333333333333
		@test I.M(2,1)==0.016666666666666663
		@test I.M(3,0)==0.05
		@test I.M(3,1)==0.008333333333333338
		@test I.M(3,2)==0.0023809523809523586
		@test I.M(3,3)==0.0008928571428571397
	end


	@testset "TT" begin
		tau=[0.0 1.0 0.0; 0.0 0.0 1.0; 0.0 0.0 0.0];
		@test I.TT(tau, 0,0,0)==0.5
		@test I.TT(tau, 1,0,0)==0.16666666666666666
		@test I.TT(tau, 1,1,0)==0.041666666666666685
		@test I.TT(tau, 1,1,1)==0.0
		@test I.TT(tau, 2,0,0)==0.08333333333333333
		@test I.TT(tau, 2,1,0)==0.016666666666666663
		@test I.TT(tau, 2,2,0)==0.005555555555555545
		@test I.TT(tau, 2,2,1)==0.0
		@test I.TT(tau, 2,2,2)==0.0
	end


	@testset "II" begin
		V = [0.0 1.0 0.0; 0.0 0.0 1.0; 0.0 0.0 0.0];
		FV = [[1,2,3]];
		P = V,FV;
		@test I.II(P, 0,0,0)==0.5
		@test I.II(P, 1,0,0)==0.16666666666666666
		@test I.II(P, 0,1,0)>=0.1666666666666666
		@test I.II(P, 0,0,1)==0.0
		@test I.II(P, 1,1,1)==0.0
	end


	@testset "III" begin
		V = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0];
		FV = [[1, 2, 4], [1, 3, 2], [4, 3, 1], [2, 3, 4]];
		P = V,FV;
		@test I.III(P, 0,0,0)>0.166666666
		@test I.III(P, 0,0,0)<0.166666888
		@test I.III(P, 1,0,0)>0.041666666
		@test I.III(P, 1,0,0)<0.041666888
		@test I.III(P, 0,1,0)>0.041666666
		@test I.III(P, 0,1,0)<0.041666888
		@test I.III(P, 0,0,1)>0.041666666
		@test I.III(P, 0,0,1)<0.041666888
		@test I.III(P, 10,10,10)>1.3377e-11
		@test I.III(P, 10,10,10)<1.3388e-11
	end


	@testset "surface" begin
		V,FV = Lar.simplexGrid([1,1]);
		P = [V;[0 0 0 0]], FV
		@test I.surface(P)==1.0
		p = Lar.Struct([Lar.t(0.5,0.5,0), Lar.r(0,0,pi/4), P]);
		q = Lar.struct2lar(p);
		q = convert(Tuple{Array{Float64,2},Array{Array{Int64,1},1}}, q)
		@test I.surface(q)>1.0000000
		@test I.surface(q)<1.0000222
	end


	@testset "volume" begin
		V = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0];
		FV = [[1, 2, 4], [1, 3, 2], [4, 3, 1], [2, 3, 4]];
		P = V,FV;
		@test I.volume(P)>0.166666666
		@test I.volume(P)<0.166668888
	end


	@testset "firstMoment" begin
		V = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0];
		FV = [[1, 2, 4], [1, 3, 2], [4, 3, 1], [2, 3, 4]];
		P = V,FV;
		@test I.firstMoment(P)[1]<0.0416667
		@test I.firstMoment(P)[1]>0.0416665

		@test I.firstMoment(P)[2]<0.0416667
		@test I.firstMoment(P)[2]>0.0416665

		@test I.firstMoment(P)[3]<0.0416667
		@test I.firstMoment(P)[3]>0.0416665

		@test abs(I.firstMoment(P)[1]-I.firstMoment(P)[2])<0.00001
		@test abs(I.firstMoment(P)[2]-I.firstMoment(P)[3])<0.00001
		@test abs(I.firstMoment(P)[3]-I.firstMoment(P)[1])<0.00001
	end


	@testset "secondMoment" begin
		V = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0];
		FV = [[1, 2, 4], [1, 3, 2], [4, 3, 1], [2, 3, 4]];
		P = V,FV;
		@test I.secondMoment(P)[1]<0.0166666669
		@test I.secondMoment(P)[1]>0.0166666664

		@test I.secondMoment(P)[2]<0.0166666669
		@test I.secondMoment(P)[2]>0.0166666664

		@test I.secondMoment(P)[3]<0.0166666669
		@test I.secondMoment(P)[3]>0.0166666664

		@test abs(I.secondMoment(P)[1]-I.secondMoment(P)[2])<0.00001
		@test abs(I.secondMoment(P)[2]-I.secondMoment(P)[3])<0.00001
		@test abs(I.secondMoment(P)[3]-I.secondMoment(P)[1])<0.00001
	end


	@testset "inertiaProduct" begin
		V = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0];
		FV = [[1, 2, 4], [1, 3, 2], [4, 3, 1], [2, 3, 4]];
		P = V,FV;
		@test I.inertiaProduct(P)[1]<0.00833666
		@test I.inertiaProduct(P)[1]>0.00833000

		@test I.inertiaProduct(P)[2]<0.00833666
		@test I.inertiaProduct(P)[2]>0.00833000

		@test I.inertiaProduct(P)[3]<0.00833666
		@test I.inertiaProduct(P)[3]>0.00833000
	end


	@testset "centroid" begin
		V = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0];
		FV = [[1, 2, 4], [1, 3, 2], [4, 3, 1], [2, 3, 4]];
		P = V,FV;
		@test I.centroid(P)[1]<0.26
		@test I.centroid(P)[1]>0.24

		@test I.centroid(P)[2]<0.26
		@test I.centroid(P)[2]>0.24

		@test I.centroid(P)[3]<0.26
		@test I.centroid(P)[3]>0.24
	end


	@testset "inertiaMoment" begin
		V = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0];
		FV = [[1, 2, 4], [1, 3, 2], [4, 3, 1], [2, 3, 4]];
		P = V,FV;
		@test I.inertiaMoment(P)[1]<0.0333555
		@test I.inertiaMoment(P)[1]>0.0333111

		@test I.inertiaMoment(P)[2]<0.0333555
		@test I.inertiaMoment(P)[2]>0.0333111

		@test I.inertiaMoment(P)[3]<0.0333555
		@test I.inertiaMoment(P)[3]>0.0333111
	end

end

