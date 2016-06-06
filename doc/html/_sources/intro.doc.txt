Introduction
=================

The objective of the Enhanced Sampling Toolkit is to provide a programming toolkit that facilitates rapid prototyping and development of enhanced sampling algorithms for use in molecular dynamics applications.

The toolkit is implemented in 100% Python and is targeted for development in the Python language. One reason for this decision is that Python is a popular language in the broader scientific community and has a wide support base for developing scientific codes. However, more importantly, a strength of the language lies in the speed at which ideas can be implemented into working code. Since rapid and flexible prototyping of new algorithms is the core priority of the toolkit, Python seems a natural choice in language. However, if the need arises in a future date, ports to other languages may be considered and integrated into the package.

In many places of the toolkit some care has been made to adhere to the powerful Python object design principles inherent in the Python data model. We find this to be a powerful feature of the Python language since it facilitates clean, "Pythonic" use of the toolkit in our applications. We are continually working to improve this aspect of the toolkit as the project develops.

In the basic sense, the toolkit serves to wrap commonly used MD codes and in the process abstract the interactions between algorithm level code and the MD engine. This abstractions provides useful extensibilty in the sense that algorithm codes that are implemented in the Walker API can be reused and swapped between MD models and even entire MD codes. Furthermore, the expensive integration steps are executed in faster compiled codes and avoids some of the inefficiencies introduced in the choice of Python.

In addition to rapid algorithm prototyping, we've found in it's development that the toolkit is effective in HPC environments as well. Because the underlying dynamics are executed in commonly used MD codes, the toolkit has access to the HPC features that have been optimized in these codes and can therefore leverage MPI parallelism concurrently at the algorithm level and the MD level as well as support for accelerators such as GPU's or Intel MIC cards.

The structure of the toolkit comprises of two parts. At the core, we provide a specification we call the Walker API. The core API serves to define a set of interactions between algorithm level code and the underlying MD engine. This specification also defines the features developers need to implement for a new dynamics engine to leverage the full features of the toolkit.

On top of this core API, we've implemented several common enhanced sampling algorithms. 