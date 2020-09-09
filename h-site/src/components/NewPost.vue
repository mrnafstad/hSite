<template>
	<div id="newpost">
		<h3>New post to the blog:</h3>
		<p><input v-model="post.title" type="string" placeholder="Title..."></p>
		<p><textarea rows="4" cols="55"
		v-model="post.postText" type="textarea" placeholder="Post..." /></p>
		<button id="subPub" @click="commitPost(post)">Submit and publish</button>
	</div>
</template>

<script>
import { mapActions } from 'vuex'
export default {
	name: "NewPost",
	data() {
		return {
			post: {
				postText: '',
				title: ''
			}
		}
	},
	methods: {
		...mapActions(['commitPost']),
	}
}
</script>

<style>
#newpost {
	border-style: outset;
	background-color: snow;
	padding: 5px;
}

#subPub {
	background-color: lightgrey;
	border-style: outset;
}
</style>