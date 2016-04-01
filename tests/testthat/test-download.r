context("Download")

test_that("Results of download shared protein network", {
    org<-SpidermiRquery_species(species)
    net_shar_prot<-SpidermiRquery_spec_networks(organismID = org[9,],
                                    network = "SHpd")
    out_net<-SpidermiRdownload_net(net_shar_prot)
	expect_that(length(out_net), equals(2))
})