function getAlgorithmParams(name::String, ϵ::Float64; newparams::Bool=true)
    if name == "numerical_expr"

        params_numerical_exprs = BBSTTDFPParams(
            0.11, # 0.98
            0.5,
            0.01,
            1.8,
            10 * ϵ,# 1e-30,
            Inf64, #1e+30,
            0.1,
            0.2,#0.5
            newparams ? 0.001 : nothing,
            newparams ? 0.6 : nothing
        )
        return params_numerical_exprs

    elseif name == "cs_expr"
        #used for cs_proccessing
        params2 = BBSTTDFPParams(
            0.9,
            0.85,
            0.9,
            1.75,
            10 * ϵ,
            Inf,
            0.94,
            0.1,
            newparams ? 0.001 : nothing,
            newparams ? 0.6 : nothing
        )

        return params2
    end
end
#used for cs_proccessing
ahdfpm_params = AHDFPMParameters(0.6, 1.75, 0.0001, 0.001, 0.4, 2.4, 0.8, 0.1, 1)

#used for cs_proccessing
cgdfp_params = CGDFPParameters(0.6, 1, 0.001, 0.7, 0.3)

#used for cs_proccessing
mopcg_params = MOPCGParameters(1, 0.6, 0.1, 0.2)

# params0 = BBSTTDFPParams(
#     0.11,
#     0.5,
#     0.01,
#     1.8,
#     1e-30,
#     1e+30,
#     0.1,
#     0.2,#0.5
#     newparams ? 0.001 : nothing,
#     newparams ? 0.6 : nothing
# )

# params1 = BBSTTDFPParams(
#     0.98,
#     0.8,
#     0.0001,
#     1.8,
#     1e-30,
#     Inf,
#     0.1,
#     0.5,
#     nothing,#0.001,
#     nothing#0.6
# )

# params1 = BBSTTDFPParams(
#     0.11408250586841473,
#     0.6570240538036174,
#     0.008575295412900363,
#     1.8928223633049506,
#     0.10912531544765902,
#     0.8312657796763663,
#     0.12207595744012179,
#     0.20929112731122723,
#     newparams ? 0.001 : nothing,
#     newparams ? 0.6 : nothing
# )
# params2 = BBSTTDFPParams(
#     0.19091055694541392,
#     0.8248500773070584,
#     0.6059149417443106,
#     1.0160210808836148,
#     0.968177563496758,
#     1.8227259069163115,
#     0.9507650425301077,
#     0.8808475746290094,
#     newparams ? 0.030826250174403214 : nothing,
#     newparams ? 0.3719729951071617 : nothing
# )
# params3 = BBSTTDFPParams(
#             0.98,
#             0.9,
#             0.0001,
#             1.8,
#             1e-30,
#             Inf,
#             0.1,
#             0.5,
#             nothing,#0.001,
#             nothing#0.6
#         )
#         params4 = BBSTTDFPParams(
#             0.19091055694541392,
#             0.8248500773070584,
#             0.6059149417443106,
#             1.0160210808836148,
#             0.968177563496758,
#             1.8227259069163115,
#             0.9507650425301077,
#             0.8808475746290094,
#             0.030826250174403214,
#             0.3719729951071617
#         )