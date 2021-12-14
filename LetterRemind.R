library(blastula)
create_smtp_creds_file(file = "haozhe_nuist",
                       user = 'haozhe_pang@nuist.edu.cn',
                       host = 'smtp.exmail.qq.com',
                       port = 465,use_ssl = TRUE)


date_time <- add_readable_time() #添加时间
email <-compose_email()
email %>%
  smtp_send(
    from = "haozhe_pang@nuist.edu.cn", # 修改发件人
    to = "haozhe_pang@nuist.edu.cn", # 收件人
    subject = "n=50已经完成",
    credentials = creds_file(file = "haozhe_nuist")
)