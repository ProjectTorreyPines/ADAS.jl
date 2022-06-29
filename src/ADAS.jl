using RedefStructs

#module ADAS


filepath = "/Users/guterl/Dropbox/api/scd93r_c.dat"
function read_adas_file(filepath::String; verbose::Bool = false)
    
    if verbose
       println("Reading ADAS file: $filepath") 
    end

    @assert isfile(filepath) "File not found: $filepath"

end
#end
@redef struct ADASHeader
    
    n_ions :: Int64
    n_ne :: Int64
    n_Te :: Int64
    imin :: Int64
    imax :: Int64
    details :: String
    function ADASHeader(header::String)
        data = parse.(Int64,split(header)[1:5])
        details =  join(split(header)[6:end]," ")
        new(data...,details)
    end
end
@redef struct ADASBlock
    header :: String
    content :: Vector{String}
end
function get_block_attr(block_header::String, attr::String)
    rg = r"$attr\s*=\s*([\d]+)"
    m = match(block_header,rg)
    return m.captures
end

function get_metastable_state()
    

function split_blocks(lines::Vector{String})
    blocks = Vector{ADASBlock}()
    iprevious = 1
    for i in 1:length(lines)
        if startswith(lines[i],"--") 
            push!(blocks, ADASBlock(lines[max(1,iprevious-1)],lines[iprevious:i-1]))
            iprevious = i+1
        end
    end
    return blocks
end


@redef struct ADASAxis
    log_Te :: Vector{Float64}
    log_ne :: Vector{Float64}
    function ADASAxis(block::Vector{String}, header :: ADASHeader)
    log_Te = Vector{Float64}()
    log_ne = Vector{Float64}()
    idx = 1
    while length(log_ne) < header.n_ne
        push!(log_ne,parse.(Float64,[n for n in split(block[idx])])...)
        idx += 1
    end
    while length(log_Te) < header.n_Te
        push!(log_Te,parse.(Float64,[t for t in split(block[idx])])...)
        idx += 1
    end
    new(log_Te, log_ne)
    end
end



@redef struct ADASRates{T<:Array{Float64,3}}
    rates :: T
    function ADASRates{T}(block::Vector{String}, header :: ADASHeader) where T
        rates = zeros(Float64,length(range(header.imin, header.imax)),header.n_ne, header.n_Te)
        for i_ion in range(header.imin, header.imax)
            blocks[i_ion]
            rates_tmp = zeros(Float64,header.n_ne * header.n_Te)
            idx = 1
            while length(rates) < header.n_ne * header.n_T
                push!(rates_tmp,parse.(Float64,[r for r in lines[idx].split()])...)
                idx += 1
            end
            rates[i_ion,:,:] = reshape(rates_tmp,n_T, n_ne)
        end
        new{T}(rates) 
    end
end
open(filepath) do fid
    # parse header
    header = readline(f)
    if verbose
        println("header : $header")
    end

    

    details = ' '.join(header.split()[3:])

    f.readline()
    n_ions, n_ne, n_T = int(n_ions), int(n_ne), int(n_T)
    imin, imax = int(imin), int(imax)
    # parse Te and ne
    logT = []
    logNe = []
    while len(logNe) < n_ne:
        line = f.readline()
        logNe = logNe + [float(n) for n in line.split()]
    while len(logT) < n_T:
        line = f.readline()
        logT = logT + [float(t) for t in line.split()]

    self.logT = np.array(logT)
    self.logNe = np.array(logNe)
    # parse data
    self.data = np.zeros((n_ions, n_T, n_ne))

    for i_ion in range(imin - 1, imax):
        f.readline()
        plsx = []
        while len(plsx) < n_ne * n_T:
            line = f.readline()
            plsx = plsx + [float(L) for L in line.split()]

        self.data[i_ion] = np.array(plsx).reshape(n_T, n_ne)
# Write your package code here.
class adas_file:
    '''Read ADAS file in ADF11 format over the given ne, T.'''

    def __init__(self, filepath):

        if not os.path.isfile(filepath):
            # backup option when is the file on localy available
            atom_data_path = '/fusion/projects/codes/strahl/public20/atomdat/newdat/'

            datafile = OMFITascii(
                atom_data_path + filepath, server=MainSettings['SERVER']['iris']['server'], tunnel=MainSettings['SERVER']['iris']['tunnel']
            )
            filepath = datafile.filename
        self.filepath = filepath
        self.filename = os.path.basename(filepath)
        self.file_type = self.filename[:3]

        if self.file_type not in ['brs', 'sxr']:
            self.imp = self.filename.split('_')[1].split('.')[0]

        # get data
        self.load()

        # settings for plotting
        self.n_ion = self.data.shape[0]
        self.ncol = np.ceil(np.sqrt(self.n_ion))
        self.nrow = np.ceil(self.n_ion / self.ncol)

    def load(self):

        with open(self.filepath) as f:
            header = f.readline()
            n_ions, n_ne, n_T, imin, imax = header.split()[:5]

            details = ' '.join(header.split()[3:])

            f.readline()
            n_ions, n_ne, n_T = int(n_ions), int(n_ne), int(n_T)
            imin, imax = int(imin), int(imax)

            logT = []
            logNe = []
            while len(logNe) < n_ne:
                line = f.readline()
                logNe = logNe + [float(n) for n in line.split()]
            while len(logT) < n_T:
                line = f.readline()
                logT = logT + [float(t) for t in line.split()]

            self.logT = np.array(logT)
            self.logNe = np.array(logNe)

            self.data = np.zeros((n_ions, n_T, n_ne))

            for i_ion in range(imin - 1, imax):
                f.readline()
                plsx = []
                while len(plsx) < n_ne * n_T:
                    line = f.readline()
                    plsx = plsx + [float(L) for L in line.split()]

                self.data[i_ion] = np.array(plsx).reshape(n_T, n_ne)

    def plot(self, fig=None, axes=None):

        if fig is None or axes is None:
            fig, axes = plt.subplots(int(self.ncol), int(self.nrow), sharex=True, sharey=True)

        axes = atleast_2d(axes)
        colormap = cm.rainbow
        fig.suptitle(self.filename + '  ' + get_file_types()[self.file_type])

        for i, ax in enumerate(axes.flatten()):
            if i >= self.n_ion:
                break
            if all(self.data[i, :, :].std(axis=1) == 0):  # independent of density
                ax.plot(self.logT, self.data[i, :, 0])
            else:
                ax.set_prop_cycle('color', colormap(np.linspace(0, 1, self.data.shape[2])))
                ax.plot(self.logT, self.data[i])
            ax.grid(True)
            if self.file_type != 'brs':
                charge = i + 1 if self.file_type in ['scd', 'prs', 'ccd', 'prb'] else i
                ax.set_title(self.imp + '$^{%d\!+}$' % (charge))

        for ax in axes[-1]:
            ax.set_xlabel('$\log T_e\ \mathrm{[eV]}$')
        for ax in axes[:, 0]:
            if self.file_type in ['scd', 'acd', 'ccd']:
                ax.set_ylabel('$\log(' + self.file_type + ')\ \mathrm{[cm^3/s]}$')
            elif self.file_type in ['prb', 'plt', 'prc', 'pls', 'brs', 'prs']:
                ax.set_ylabel('$\log(' + self.file_type + ')\ \mathrm{[W\cdot cm^3]}$')
end
